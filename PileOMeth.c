#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>

inline int isCpG(char *seq, int pos, int seqlen) {
    if(pos+1 == seqlen) return 0;
    if(*(seq+pos) == 'C' || *(seq+pos) == 'c') {
        if(*(seq+pos+1) == 'G' || *(seq+pos+1) == 'G') return 1;
        return 0;
    } else if(*(seq+pos) == 'G' || *(seq+pos) == 'g') {
        if(pos == 0) return 0;
        if(*(seq+pos-1) == 'C' || *(seq+pos-1) == 'c') return -1;
        return 0;
    }
    return 0;
}

inline int isCHG(char *seq, int pos, int seqlen) {
    if(pos+2 >= seqlen) return 0;
    if(*(seq+pos) == 'C' || *(seq+pos) == 'c') {
        if(*(seq+pos+2) == 'G' || *(seq+pos+2) == 'G') return 1;
        return 0;
    } else if(*(seq+pos) == 'G' || *(seq+pos) == 'g') {
        if(pos <= 1) return 0;
        if(*(seq+pos-2) == 'C' || *(seq+pos-2) == 'c') return -1;
        return 0;
    }
    return 0;
}

inline int isCHH(char *seq, int pos, int seqlen) {
    if(pos+2 >= seqlen) return 0;
    if(*(seq+pos) == 'C' || *(seq+pos) == 'c') return 1;
    else if(*(seq+pos) == 'G' || *(seq+pos) == 'g') {
        if(pos <= 1) return 0;
        return -1;
    }
    return 0;
}

typedef struct {
    int KeepCpG, KeepCHG, KeepCHH;
    int minMapq, minPhred, ignoreDupes, maxDepth;
    FILE **output_fp;
    char *reg;
    htsFile *fp;
    hts_idx_t *bai;
    faidx_t *fai;
} Config;

typedef struct {
    Config *config;
    bam_hdr_t *hdr;
    hts_itr_t *iter;
} mplp_data;

int getStrand(bam1_t *b) {
    if(b->core.flag & BAM_FPAIRED) {
        if((b->core.flag & 0x50) == 0x50) return 1; //Read1, reverse comp. == OB
        else if(b->core.flag & 0x40) return 0; //Read1, forward == OT
        else if((b->core.flag & 0x90) == 0x90) return 0; //Read2, reverse comp. == OT
        else if(b->core.flag & 0x80) return 0; //Read2, forward == OB
        return -1; //One of the above should be set!
    } else {
        if(b->core.flag & 0x10) return 1; //Reverse comp. == OB
        return 0;
    }
}

int updateMetrics(Config *config, const bam_pileup1_t *plp) {
    uint8_t base = bam_seqi(bam_get_seq(plp->b), plp->qpos);
    int strand = getStrand(plp->b); //0=OT, 1=OB, 2=CTOT, 3=CTOB

    //Is the phred score even high enough?
    if(bam_get_qual(plp->b)[plp->qpos] < config->minMapq) return 0;

    if(base == 2 && strand == 0) return 1; //C on an OT read
    else if(base == 8 && strand == 0)  return -1; //T on an OT read
    else if(base == 4 && strand == 1) return 1; //G on an OB read
    else if(base == 1 && strand == 1) return -1; //A on an OB read
    return 0;
}

//This needs to be defined
int filter_func(void *data, bam1_t *b) {
    int rv, NH;
    mplp_data *ldata = (mplp_data *) data;
    uint8_t *p;

    while(1) {
        rv = ldata->iter ? sam_itr_next(ldata->config->fp, ldata->iter, b) : sam_read1(ldata->config->fp, ldata->hdr, b);
        if(rv<0) return rv;
        if(b->core.tid == -1 || b->core.flag & BAM_FUNMAP) continue;
        if(b->core.qual < ldata->config->minMapq) continue;
        if(b->core.flag & 0x300) continue; //Ignore secondary alignments and those with QC failed
        if(ldata->config->ignoreDupes && b->core.flag & BAM_FDUP) continue;
        p = bam_aux_get(b, "NH");
        if(p != NULL) {
            NH = bam_aux2i(p);
            if(NH>1) continue; //Ignore obvious multimappers
        }
        break;
    }
    return rv;
}

void extractCalls(Config *config) {
    bam_hdr_t *hdr = sam_hdr_read(config->fp);
    bam_mplp_t iter;
    int ret, tid, pos, i, seqlen, type, rv;
    int beg0 = 0, end0 = 1u<<29;
    int n_plp; //This will need to be modified for multiple input files
    int ctid = -1; //The tid of the contig whose sequence is stored in "seq"
    uint32_t nmethyl, nunmethyl;
    const bam_pileup1_t **plp = calloc(1, sizeof(bam_pileup1_t *)); //This will have to be modified for multiple input files
    char *seq = NULL;
    mplp_data *data = malloc(sizeof(mplp_data));

    data->config = config;
    data->hdr = hdr;
    if (config->reg) {
        if((data->iter = sam_itr_querys(config->bai, hdr, config->reg)) == 0){
            fprintf(stderr, "failed to parse regions %s", config->reg);
            exit(1);
        }
   }

    //Start the pileup
    iter = bam_mplp_init(1, filter_func, (void **) &data);
    bam_mplp_init_overlaps(iter);
    bam_mplp_set_maxcnt(iter, config->maxDepth);
    while((ret = bam_mplp_auto(iter, &tid, &pos, &n_plp, plp)) > 0) {
        //Do we need to process this position?
	if (config->reg && (pos < beg0 || pos >= end0)) continue; // out of the region requested
        if(tid != ctid) {
            if(seq != NULL) free(seq);
            seq = faidx_fetch_seq(config->fai, hdr->target_name[tid], 0, faidx_seq_len(config->fai, hdr->target_name[tid]), &seqlen);
            if(seqlen < 0) {
                fprintf(stderr, "faidx_fetch_seq returned %i while trying to fetch the sequence for tid %i (%s)!\n",\
                    seqlen, tid, hdr->target_name[tid]);
                fprintf(stderr, "Note that the output will be truncated!\n");
                return;
            }
            ctid = tid;
        }
        if(config->KeepCpG && isCpG(seq, pos, seqlen)) {
            type = 0;
        } else if(config->KeepCHG && isCHG(seq, pos, seqlen)) {
            type = 1;
        } else if(config->KeepCHH && isCHH(seq, pos, seqlen)) {
            type = 2;
        } else {
            continue;
        }

        nmethyl = nunmethyl = 0;
        for(i=0; i<n_plp; i++) {
            rv = updateMetrics(config, plp[0]+i);
            if(rv > 0) nmethyl++;
            else if(rv<0) nunmethyl++;
        }

        if(nmethyl+nunmethyl) fprintf(config->output_fp[type], "%s\t%i\t%i\t%f\t%" PRIu32 "\t%" PRIu32 "\n", \
            hdr->target_name[tid], pos, pos+1, 1000.0 * ((double) nmethyl)/(nmethyl+nunmethyl), nmethyl, nunmethyl);
    }

    bam_hdr_destroy(hdr);
    bam_mplp_destroy(iter);
    free(data);
    free(plp);
    if(seq != NULL) free(seq);
}

void usage(char *prog) {
    fprintf(stderr, "\nUsage: %s [OPTIONS] reference.fa file.bam\n", prog);
    fprintf(stderr,
"\n"
"Options:\n"
" -q INT           Minimum MAPQ threshold to include an alignment (default 5)\n"
" -p INT           Minimum Phred threshold to include a base (default 10)\n"
" -D INT           Maximum per-base depth (default 2000)\n"
" -r STR           Region string in which to extract methylation\n"
" -o, --opref STR  Output filename prefix. CpG/CHG/CHH metrics will be\n"
"                  output to STR_CpG.bedGraph and so on.\n"
" --KeepDupes      By default, any alignment marked as a duplicate is ignored.\n"
"                  This option causes them to be incorporated.\n"
" --noCpG          Do not output CpG methylation metrics\n"
" --CHG            Output CHG methylation metrics\n"
" --CHH            Output CHH methylation metrics\n");
}

int main(int argc, char *argv[]) {
    char *opref = NULL, *oname, *p;
    int c;
    Config config;

    //Defaults
    config.KeepCpG = 1; config.KeepCHG = 0; config.KeepCHH = 0;
    config.minMapq = 5; config.minPhred = 5; config.ignoreDupes = 1;
    config.maxDepth = 2000;
    config.fai = NULL;
    config.fp = NULL;
    config.bai = NULL;
    config.reg = NULL;

    static struct option lopts[] = {
        {"opref", 1, NULL,              'o'},
        {"noCpG", 0, NULL,                1},
        {"CHG",   0, NULL,                2},
        {"CHH",   0, NULL,                3},
        {"ignoreDupes", 0, NULL,              '4'},
        {"region",      1, NULL,              'r'},
        {"help",  0, NULL,              'h'}
    };
    while((c = getopt_long(argc, argv, "q:p:r:o:D:", lopts,NULL)) >= 0) {
        switch(c) {
        case 'h' :
            usage(argv[0]);
            return 0;
        case 'o' :
            opref = strdup(optarg);
            break;
        case 'D' :
            config.maxDepth = atoi(optarg);
            break;
	case 'r':
	    config.reg = strdup(optarg);
	    break;
        case 1 :
            config.KeepCpG = 0;
            break;
        case 2 :
            config.KeepCHG = 1;
            break;
        case 3 :
            config.KeepCHH = 1;
            break;
        case 4 :
            config.ignoreDupes = 0;
            break;
        case 'q' :
            config.minMapq = atoi(optarg);
            break;
        case 'p' :
            config.minPhred = atoi(optarg);
            break;
        default :
            fprintf(stderr, "Invalid option '%c'\n", c);
            usage(argv[0]);
            return 1;
        }
    }

    if(argc == 1) {
        usage(argv[0]);
        return 0;
    }
    if(argc-optind != 2) {
        fprintf(stderr, "You must supply a reference genome in fasta format and an input BAM file!!!\n");
        usage(argv[0]);
        return -1;
    }

    //Open the files
    config.fai = fai_load(argv[optind]);
    if(config.fai == NULL) {
        fprintf(stderr, "Couldn't open the index for %s!\n", argv[optind]);
        usage(argv[0]);
        return -2;
    }
    config.fp = hts_open(argv[optind+1], "rb"); //Does this return NULL on error?
    if(config.fp == NULL) {
        fprintf(stderr, "Couldn't open %s for reading!\n", argv[optind+1]);
        return -4;
    }
    config.bai = sam_index_load(config.fp, argv[optind+1]);
    if(config.bai == NULL) {
        fprintf(stderr, "Couldn't open the index for %s!\n", argv[optind+1]);
        return -5;
    }

    //Output files
    config.output_fp = malloc(sizeof(FILE *) * 3);
    if(opref == NULL) {
        opref = strdup(argv[optind+1]);
        p = strrchr(opref, '.');
        *p = '\0';
        fprintf(stderr, "writing to prefix:'%s'\n", opref);
    }
    oname = malloc(sizeof(char) * (strlen(opref)+14));
    if(config.KeepCpG) {
        sprintf(oname, "%s_CpG.bedGraph", opref);
        config.output_fp[0] = fopen(oname, "w");
        if(config.output_fp[0] == NULL) {
            fprintf(stderr, "Couldn't open the output CpG metrics file for writing! Insufficient permissions?\n");
            return -3;
        }
        fprintf(config.output_fp[0], "track type=\"bedGraph\" description=\"%s CpG methylation levels\"\n", opref);
    }
    if(config.KeepCHG) {
        sprintf(oname, "%s_CHG.bedGraph", opref);
        config.output_fp[1] = fopen(oname, "w");
        if(config.output_fp[1] == NULL) {
            fprintf(stderr, "Couldn't open the output CHG metrics file for writing! Insufficient permissions?\n");
            return -3;
        }
        fprintf(config.output_fp[0], "track type=\"bedGraph\" description=\"%s CHG methylation levels\"\n", opref);
    }
    if(config.KeepCHH) {
        sprintf(oname, "%s_CHH.bedGraph", opref);
        config.output_fp[2] = fopen(oname, "w");
        if(config.output_fp[2] == NULL) {
            fprintf(stderr, "Couldn't open the output CHH metrics file for writing! Insufficient permissions?\n");
            return -3;
        }
        fprintf(config.output_fp[0], "track type=\"bedGraph\" description=\"%s CHH methylation levels\"\n", opref);
    }

    //Run the pileup
    extractCalls(&config);

    //Close things up
    hts_close(config.fp);
    fai_destroy(config.fai);
    if(config.KeepCpG) fclose(config.output_fp[0]);
    if(config.KeepCHG) fclose(config.output_fp[1]);
    if(config.KeepCHH) fclose(config.output_fp[2]);
    hts_idx_destroy(config.bai);
    free(opref);
    free(oname);
    free(config.output_fp);

    return 0;
}
