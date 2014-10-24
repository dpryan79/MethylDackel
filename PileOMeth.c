#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>
#include "PileOMeth.h"

inline int isCpG(char *seq, int pos, int seqlen) {
    if(pos+1 == seqlen) return 0;
    if(*(seq+pos) == 'C' || *(seq+pos) == 'c') {
        if(*(seq+pos+1) == 'G' || *(seq+pos+1) == 'g') return 1;
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
        if(*(seq+pos+2) == 'G' || *(seq+pos+2) == 'g') return 1;
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
    Config *config;
    bam_hdr_t *hdr;
    hts_itr_t *iter;
} mplp_data;

int getStrand(bam1_t *b) {
    char *XG = (char *) bam_aux_get(b, "XG");
    if(XG == NULL) { //Can't handle non-directional libraries!
        if(b->core.flag & BAM_FPAIRED) {
            if((b->core.flag & 0x50) == 0x50) return 2; //Read1, reverse comp. == OB
            else if(b->core.flag & 0x40) return 1; //Read1, forward == OT
            else if((b->core.flag & 0x90) == 0x90) return 1; //Read2, reverse comp. == OT
            else if(b->core.flag & 0x80) return 2; //Read2, forward == OB
            return 0; //One of the above should be set!
        } else {
            if(b->core.flag & 0x10) return 2; //Reverse comp. == OB
            return 0;
        }
    } else {
        if(*(XG+1) == 'C') { //OT or CTOT, due to C->T converted genome
            if((b->core.flag & 0x51) == 0x41) return 1; //Read#1 forward == OT
            else if((b->core.flag & 0x51) == 0x51) return 3; //Read #1 reverse == CTOT
            else if((b->core.flag & 0x91) == 0x81) return 3; //Read #2 forward == CTOT
            else if((b->core.flag & 0x91) == 0x91) return 1; //Read #2 reverse == OT
            else if(b->core.flag & 0x10) return 3; //Single-end reverse == CTOT
            else return 1; //Single-end forward == OT
        } else {
            if((b->core.flag & 0x51) == 0x41) return 4; //Read#1 forward == CTOB
            else if((b->core.flag & 0x51) == 0x51) return 2; //Read #1 reverse == OB
            else if((b->core.flag & 0x91) == 0x81) return 2; //Read #2 forward == OB
            else if((b->core.flag & 0x91) == 0x91) return 4; //Read #2 reverse == CTOB
            else if(b->core.flag & 0x10) return 2; //Single-end reverse == OB
            else return 4; //Single-end forward == CTOB
        }
    }
}

int updateMetrics(Config *config, const bam_pileup1_t *plp) {
    uint8_t base = bam_seqi(bam_get_seq(plp->b), plp->qpos);
    int strand = getStrand(plp->b); //1=OT, 2=OB, 3=CTOT, 4=CTOB

    assert(("Can't determine the strand of a read!", strand != 0));

    //Is the phred score even high enough?
    if(bam_get_qual(plp->b)[plp->qpos] < config->minPhred) return 0;

    if(base == 2 && (strand & 1)) return 1; //C on an OT/CTOT alignment
    else if(base == 8 && (strand & 1))  return -1; //T on an OT/CTOT alignment
    else if(base == 4 && (strand & 2)) return 1; //G on an OB/CTOB alignment
    else if(base == 1 && (strand & 2)) return -1; //A on an OB/CTOB alignment
    return 0;
}

//This will need to be restructured to handle multiple input files
int filter_func(void *data, bam1_t *b) {
    int rv, NH, overlap;
    static int32_t idxBED = 0;
    mplp_data *ldata = (mplp_data *) data;
    uint8_t *p;

    while(1) {
        rv = ldata->iter ? sam_itr_next(ldata->config->fp, ldata->iter, b) : sam_read1(ldata->config->fp, ldata->hdr, b);

        if(rv<0) return rv;
        if(b->core.tid == -1 || b->core.flag & BAM_FUNMAP) continue; //Unmapped
        if(b->core.qual < ldata->config->minMapq) continue; //-q
        if(b->core.flag & 0x300) continue; //Ignore secondary alignments and those with QC failed
        if(!ldata->config->keepDupes && b->core.flag & BAM_FDUP) continue;
        p = bam_aux_get(b, "NH");
        if(p != NULL) {
            NH = bam_aux2i(p);
            if(NH>1) continue; //Ignore obvious multimappers
        }
        if(!ldata->config->keepSingleton && (b->core.flag & 0x9) == 0x9) continue; //Singleton
        if(!ldata->config->keepDiscordant && (b->core.flag & 0x3) == 0x1) continue; //Discordant
        if((b->core.flag & 0x9) == 0x1) b->core.flag |= 0x2; //Discordant pairs can cause double counts
        if(ldata->config->bed) { //Prefilter reads overlapping a BED file (N.B., strand independent).
            overlap = spanOverlapsBED(b->core.tid, b->core.pos, bam_endpos(b), ldata->config->bed, &idxBED);
            if(overlap == 0) continue;
            if(overlap < 0) {
                rv = -1;
                break;
            }
        }
        break;
    }
    return rv;
}

void extractCalls(Config *config) {
    bam_hdr_t *hdr = sam_hdr_read(config->fp);
    bam_mplp_t iter;
    int ret, tid, pos, i, seqlen, type, rv, o = 0;
    int beg0 = 0, end0 = 1u<<29;
    int n_plp; //This will need to be modified for multiple input files
    int ctid = -1; //The tid of the contig whose sequence is stored in "seq"
    int idxBED = 0;
    uint32_t nmethyl, nunmethyl;
    const bam_pileup1_t **plp = NULL;
    char *seq = NULL, base;
    mplp_data *data = NULL;

    data = calloc(1,sizeof(mplp_data));
    if (config->reg) {
        if((data->iter = sam_itr_querys(config->bai, hdr, config->reg)) == 0){
            fprintf(stderr, "failed to parse regions %s", config->reg);
            return;
        }
    }
    if(config->bedName) {
        config->bed = parseBED(config->bedName, hdr);
        if(config->bed == NULL) return;
    }

    if(data == NULL) {
        fprintf(stderr, "Couldn't allocate space for the data structure in extractCalls()!\n");
        return;
    }
    data->config = config;
    data->hdr = hdr;
    plp = calloc(1, sizeof(bam_pileup1_t *)); //This will have to be modified for multiple input files
    if(plp == NULL) {
        fprintf(stderr, "Couldn't allocate space for the plp structure in extractCalls()!\n");
        return;
    }

    //Start the pileup
    iter = bam_mplp_init(1, filter_func, (void **) &data);
    bam_mplp_init_overlaps(iter);
    bam_mplp_set_maxcnt(iter, config->maxDepth);
    while((ret = bam_mplp_auto(iter, &tid, &pos, &n_plp, plp)) > 0) {
        //Do we need to process this position?
	if (config->reg){
	    beg0 = data->iter->beg, end0 = data->iter->end;
	    if ((pos < beg0 || pos >= end0)) continue; // out of the region requested
	}
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

        if(config->bed) { //Handle -l
            while((o = posOverlapsBED(tid, pos, config->bed, idxBED)) == -1) idxBED++;
            if(o == 0) continue; //Wrong strand
        }

        if(isCpG(seq, pos, seqlen)) {
            if(!config->keepCpG) continue;
            type = 0;
        } else if(isCHG(seq, pos, seqlen)) {
            if(!config->keepCHG) continue;
            type = 1;
        } else if(isCHH(seq, pos, seqlen)) {
            if(!config->keepCHH) continue;
            type = 2;
        } else {
            continue;
        }

        nmethyl = nunmethyl = 0;
        base = *(seq+pos);
        for(i=0; i<n_plp; i++) {
            if(config->bed) if(!readStrandOverlapsBED(plp[0][i].b, config->bed->region[idxBED])) continue;
            if(getStrand((plp[0]+i)->b) & 1) {
                if(base != 'C' && base != 'c') continue;
            } else {
                if(base != 'G' && base != 'g') continue;
            }
            rv = updateMetrics(config, plp[0]+i);
            if(rv > 0) nmethyl++;
            else if(rv<0) nunmethyl++;
        }

        if(nmethyl+nunmethyl) fprintf(config->output_fp[type], "%s\t%i\t%i\t%i\t%" PRIu32 "\t%" PRIu32 "\n", \
            hdr->target_name[tid], pos, pos+1, (int) (1000.0 * ((double) nmethyl)/(nmethyl+nunmethyl)), nmethyl, nunmethyl);
    }

    bam_hdr_destroy(hdr);
    if(data->iter) hts_itr_destroy(data->iter);
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
" -p INT           Minimum Phred threshold to include a base (default 10). This\n"
"                  must be >0.\n"
" -D INT           Maximum per-base depth (default 2000)\n"
" -r STR           Region string in which to extract methylation\n"
" -l FILE          A BED file listing regions for inclusion. Note that unlike\n"
"                  samtools mpileup, this option will utilize the strand column\n"
"                  (column 6) if present. Thus, if a region has a '+' in this\n"
"                  column, then only metrics from the top strand will be\n"
"                  output. Note that the -r option can be used to limit the\n"
"                  regions of -l.\n"
" -o, --opref STR  Output filename prefix. CpG/CHG/CHH metrics will be\n"
"                  output to STR_CpG.bedGraph and so on.\n"
" --keepDupes      By default, any alignment marked as a duplicate is ignored.\n"
"                  This option causes them to be incorporated.\n"
" --keepSingleton  By default, if only one read in a pair aligns (a singleton)\n"
"                  then it's ignored.\n"
" --keepDiscordant By default, paired-end alignments with the properly-paired bit\n"
"                  unset in the FLAG field are ignored. Note that the definition\n"
"                  of concordant and discordant is based on your aligner\n"
"                  settings.\n"
" --noCpG          Do not output CpG methylation metrics\n"
" --CHG            Output CHG methylation metrics\n"
" --CHH            Output CHH methylation metrics\n");
}

int main(int argc, char *argv[]) {
    char *opref = NULL, *oname, *p;
    int c;
    Config config;

    //Defaults
    config.keepCpG = 1; config.keepCHG = 0; config.keepCHH = 0;
    config.minMapq = 5; config.minPhred = 10; config.keepDupes = 0;
    config.keepSingleton = 0, config.keepDiscordant = 0;
    config.maxDepth = 2000;
    config.fai = NULL;
    config.fp = NULL;
    config.bai = NULL;
    config.reg = NULL;
    config.bedName = NULL;
    config.bed = NULL;

    static struct option lopts[] = {
        {"opref",        1, NULL, 'o'},
        {"noCpG",        0, NULL,   1},
        {"CHG",          0, NULL,   2},
        {"CHH",          0, NULL,   3},
        {"keepDupes",    0, NULL,   4},
        {"keepSingleton",  0, NULL,   5},
        {"keepDiscordant", 0, NULL,   6},
        {"help",         0, NULL, 'h'}
    };
    while((c = getopt_long(argc, argv, "q:p:r:l:o:D:", lopts,NULL)) >= 0) {
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
        case 'l' :
            config.bedName = optarg;
            break;
        case 1 :
            config.keepCpG = 0;
            break;
        case 2 :
            config.keepCHG = 1;
            break;
        case 3 :
            config.keepCHH = 1;
            break;
        case 4 :
            config.keepDupes = 1;
            break;
        case 5 :
            config.keepSingleton = 1;
            break;
        case 6 :
            config.keepDiscordant = 1;
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

    //Are the options reasonable?
    if(config.minPhred < 1) {
        fprintf(stderr, "-p %i is invalid. resetting to 1, which is the lowest possible value.\n", config.minPhred);
        config.minPhred = 1;
    }
    if(config.minMapq < 0) {
        fprintf(stderr, "-q %i is invalid. Resetting to 0, which is the lowest possible value.\n", config.minMapq);
        config.minMapq = 0;
    }

    //Is there still a metric to output?
    if(!(config.keepCpG + config.keepCHG + config.keepCHH)) {
        fprintf(stderr, "You haven't specified any metrics to output!\nEither don't use the --noCpG option or specify --CHG and/or --CHH.\n");
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
        if(p != NULL) *p = '\0';
        fprintf(stderr, "writing to prefix:'%s'\n", opref);
    }
    oname = malloc(sizeof(char) * (strlen(opref)+14));
    if(config.keepCpG) {
        sprintf(oname, "%s_CpG.bedGraph", opref);
        config.output_fp[0] = fopen(oname, "w");
        if(config.output_fp[0] == NULL) {
            fprintf(stderr, "Couldn't open the output CpG metrics file for writing! Insufficient permissions?\n");
            return -3;
        }
        fprintf(config.output_fp[0], "track type=\"bedGraph\" description=\"%s CpG methylation levels\"\n", opref);
    }
    if(config.keepCHG) {
        sprintf(oname, "%s_CHG.bedGraph", opref);
        config.output_fp[1] = fopen(oname, "w");
        if(config.output_fp[1] == NULL) {
            fprintf(stderr, "Couldn't open the output CHG metrics file for writing! Insufficient permissions?\n");
            return -3;
        }
        fprintf(config.output_fp[1], "track type=\"bedGraph\" description=\"%s CHG methylation levels\"\n", opref);
    }
    if(config.keepCHH) {
        sprintf(oname, "%s_CHH.bedGraph", opref);
        config.output_fp[2] = fopen(oname, "w");
        if(config.output_fp[2] == NULL) {
            fprintf(stderr, "Couldn't open the output CHH metrics file for writing! Insufficient permissions?\n");
            return -3;
        }
        fprintf(config.output_fp[2], "track type=\"bedGraph\" description=\"%s CHH methylation levels\"\n", opref);
    }

    //Run the pileup
    extractCalls(&config);

    //Close things up
    hts_close(config.fp);
    fai_destroy(config.fai);
    if(config.keepCpG) fclose(config.output_fp[0]);
    if(config.keepCHG) fclose(config.output_fp[1]);
    if(config.keepCHH) fclose(config.output_fp[2]);
    hts_idx_destroy(config.bai);
    free(opref);
    if(config.reg) free(config.reg);
    if(config.bed) destroyBED(config.bed);
    free(oname);
    free(config.output_fp);

    return 0;
}
