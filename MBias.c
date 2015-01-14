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

//Returns 1 if methylated, -1 if unmethylated and 0 otherwise
int updateMetrics(Config *config, const bam_pileup1_t *plp) {
    uint8_t base = bam_seqi(bam_get_seq(plp->b), plp->qpos);
    int strand = getStrand(plp->b); //1=OT, 2=OB, 3=CTOT, 4=CTOB

    if(strand==0) {
        fprintf(stderr, "Can't determine the strand of a read!\n");
        assert(strand != 0);
    }

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
    int rv, NH;
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
        break;
    }
    return rv;
}

strandMeth *growStrandMeth(strandMeth *s, int32_t l) {
    int32_t m = kroundup32(l);
    int i;
    if(m==0) m=32; //Enforce a minimum length

    s->unmeth1 = realloc(s->unmeth1, sizeof(uint32_t)*m);
    s->meth1 = realloc(s->meth1, sizeof(uint32_t)*m);
    s->unmeth2 = realloc(s->unmeth2, sizeof(uint32_t)*m);
    s->meth2 = realloc(s->meth2, sizeof(uint32_t)*m);
    assert(s->unmeth1);
    assert(s->unmeth2);
    assert(s->meth1);
    assert(s->meth2);

    for(i=s->m; i<m; i++) {
        s->unmeth1[i] = 0;
        s->meth1[i] = 0;
        s->unmeth2[i] = 0;
        s->meth2[i] = 0;
    }
    s->m = m;
    return s;
}

void extractCalls(Config *config, char *opref, int SVG, int txt) {
    bam_hdr_t *hdr = sam_hdr_read(config->fp);
    bam_mplp_t iter;
    int ret, tid, pos, i, seqlen, rv;
    int strand;
    int n_plp; //This will need to be modified for multiple input files
    int ctid = -1; //The tid of the contig whose sequence is stored in "seq"
    const bam_pileup1_t **plp = NULL;
    char *seq = NULL, base;
    mplp_data *data = NULL;
    strandMeth *meths[4];

    for(i=0; i<4; i++) {
        meths[i] = calloc(1, sizeof(strandMeth));
        assert(meths[i]);
    }

    data = calloc(1,sizeof(mplp_data));
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

        if(isCpG(seq, pos, seqlen)) {
            if(!config->keepCpG) continue;
        } else if(isCHG(seq, pos, seqlen)) {
            if(!config->keepCHG) continue;
        } else if(isCHH(seq, pos, seqlen)) {
            if(!config->keepCHH) continue;
        } else {
            continue;
        }

        base = *(seq+pos);
        for(i=0; i<n_plp; i++) {
            strand = getStrand((plp[0]+i)->b);
            if(strand & 1) {
                if(base != 'C' && base != 'c') continue;
            } else {
                if(base != 'G' && base != 'g') continue;
            }
            rv = updateMetrics(config, plp[0]+i);
            if(rv != 0) {
                if((plp[0]+i)->qpos >= meths[strand-1]->m)
                    meths[strand-1] = growStrandMeth(meths[strand-1], (plp[0]+i)->qpos);
                if(rv > 0) {
                    if((plp[0]+i)->b->core.flag & BAM_FREAD2) {
                        assert((meths[strand-1]->unmeth2[(plp[0]+i)->qpos]) < 0xFFFFFFFF);
                        meths[strand-1]->unmeth2[(plp[0]+i)->qpos]++;
                    } else {
                        assert((meths[strand-1]->unmeth1[(plp[0]+i)->qpos]) < 0xFFFFFFFF);
                        meths[strand-1]->unmeth1[(plp[0]+i)->qpos]++;
                    }
                } else {
                    if((plp[0]+i)->b->core.flag & BAM_FREAD2) {
                        assert((meths[strand-1]->meth2[(plp[0]+i)->qpos]) < 0xFFFFFFFF);
                        meths[strand-1]->meth2[(plp[0]+i)->qpos]++;
                    } else {
                        assert((meths[strand-1]->meth1[(plp[0]+i)->qpos]) < 0xFFFFFFFF);
                        meths[strand-1]->meth1[(plp[0]+i)->qpos]++;
                    }
                }
                if((plp[0]+i)->qpos > meths[strand-1]->l) meths[strand-1]->l = (plp[0]+i)->qpos;
            }
        }
    }

    //Report some output
    if(SVG) makeSVGs(opref, meths);
    if(txt) makeTXT(meths);

    //Clean up
    bam_hdr_destroy(hdr);
    if(data->iter) hts_itr_destroy(data->iter);
    bam_mplp_destroy(iter);
    free(data);
    free(plp);
    if(seq != NULL) free(seq);
    for(i=0; i<4; i++) {
        if(meths[i]->meth1) free(meths[i]->meth1);
        if(meths[i]->unmeth1) free(meths[i]->unmeth1);
        if(meths[i]->meth2) free(meths[i]->meth2);
        if(meths[i]->unmeth2) free(meths[i]->unmeth2);
        free(meths[i]);
    }
}

void usage(char *prog) {
    fprintf(stderr, "\nUsage: %s [OPTIONS] reference.fa file.bam prefix\n", prog);
    fprintf(stderr,
"\n"
"Options:\n"
" -q INT           Minimum MAPQ threshold to include an alignment (default 5)\n"
" -p INT           Minimum Phred threshold to include a base (default 10). This\n"
"                  must be >0.\n"
" -D INT Maximum per-base depth (default 2000)\n"
" --keepDupes      By default, any alignment marked as a duplicate is ignored.\n"
"                  This option causes them to be incorporated.\n"
" --keepSingleton  By default, if only one read in a pair aligns (a singleton)\n"
"                  then it's ignored.\n"
" --keepDiscordant By default, paired-end alignments with the properly-paired bit\n"
"                  unset in the FLAG field are ignored. Note that the definition\n"
"                  of concordant and discordant is based on your aligner\n"
"                  settings.\n"
" --txt            Output tab separated metrics to the screen. These can be\n"
"                  imported into R or another program for manual plotting and\n"
"                  analysis.\n"
" --noSVG          Don't produce the SVG files. This option implies --txt. Note\n"
"                  that an output prefix is no longer required with this option.\n"
" --noCpG          Do not output CpG methylation metrics\n"
" --CHG            Output CHG methylation metrics\n"
" --CHH            Output CHH methylation metrics\n");
}

int main(int argc, char *argv[]) {
    char *opref = NULL;
    int c, SVG = 1, txt = 0;
    Config config;

    //Defaults
    config.keepCpG = 1; config.keepCHG = 0; config.keepCHH = 0;
    config.minMapq = 5; config.minPhred = 10; config.keepDupes = 0;
    config.keepSingleton = 0, config.keepDiscordant = 0;
    config.fai = NULL;
    config.fp = NULL;
    config.bai = NULL;
    config.reg = NULL;
    config.bedName = NULL;
    config.bed = NULL;

    static struct option lopts[] = {
        {"noCpG",        0, NULL,   1},
        {"CHG",          0, NULL,   2},
        {"CHH",          0, NULL,   3},
        {"keepDupes",    0, NULL,   4},
        {"keepSingleton",  0, NULL,   5},
        {"keepDiscordant", 0, NULL,   6},
        {"txt",          0, NULL,   7},
        {"noSVG",        0, NULL,   8},
        {"help",         0, NULL, 'h'}
    };
    while((c = getopt_long(argc, argv, "q:p:D:", lopts,NULL)) >= 0) {
        switch(c) {
        case 'h' :
            usage(argv[0]);
            return 0;
        case 'D' :
            config.maxDepth = atoi(optarg);
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
        case 7 :
            txt = 1;
            break;
        case 8 :
            SVG = 0;
            txt = 1;
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
    if((SVG && argc-optind != 3) || (!SVG && argc-optind < 2)) {
        fprintf(stderr, "You must supply a reference genome in fasta format, an input BAM file, and an output prefix!!!\n");
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
    if((config.fai = fai_load(argv[optind])) == NULL) {
        fprintf(stderr, "Couldn't open the index for %s!\n", argv[optind]);
        usage(argv[0]);
        return -2;
    }
    if((config.fp = hts_open(argv[optind+1], "rb")) == NULL) {
        fprintf(stderr, "Couldn't open %s for reading!\n", argv[optind+1]);
        return -4;
    }
    if((config.bai = sam_index_load(config.fp, argv[optind+1])) == NULL) {
        fprintf(stderr, "Couldn't load the index for %s, will attempt to build it.\n", argv[optind+1]);
        if(bam_index_build(argv[optind+1], 0) < 0) {
            fprintf(stderr, "Couldn't build the index for %s! File corrupted?\n", argv[optind+1]);
            return -5;
        }
        if((config.bai = sam_index_load(config.fp, argv[optind+1])) == NULL) {
            fprintf(stderr, "Still couldn't load the index, quiting.\n");
            return -5;
        }
    }

    //Output files (this needs to be filled in)
    if(SVG) opref = argv[optind+2];

    //Run the pileup
    extractCalls(&config, opref, SVG, txt);

    //Close things up
    hts_close(config.fp);
    fai_destroy(config.fai);
    hts_idx_destroy(config.bai);

    return 0;
}
