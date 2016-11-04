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

void print_version(void);

strandMeth *growStrandMeth(strandMeth *s, int32_t l) {
    int32_t m;
    int i;
    l++;
    m = kroundup32(l);
    if(m<32) m=32; //Enforce a minimum length

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

void extractMBias(Config *config, char *opref, int SVG, int txt) {
    bam_hdr_t *hdr = sam_hdr_read(config->fp);
    bam_mplp_t iter;
    int ret, tid, pos, i, seqlen, rv, o = 0;
    int beg0 = 0, end0 = 1u<<29;
    int strand;
    int n_plp; //This will need to be modified for multiple input files
    int ctid = -1; //The tid of the contig whose sequence is stored in "seq"
    int idxBED = 0;
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


    plp = calloc(1, sizeof(bam_pileup1_t *)); //This will have to be modified for multiple input files
    if(plp == NULL) {
        fprintf(stderr, "Couldn't allocate space for the plp structure in extractCalls()!\n");
        return;
    }

    //Start the pileup
    iter = bam_mplp_init(1, filter_func, (void **) &data);
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
        } else if(isCHG(seq, pos, seqlen)) {
            if(!config->keepCHG) continue;
        } else if(isCHH(seq, pos, seqlen)) {
            if(!config->keepCHH) continue;
        } else {
            continue;
        }

        base = *(seq+pos);
        for(i=0; i<n_plp; i++) {
            if(config->bed) if(!readStrandOverlapsBED(plp[0][i].b, config->bed->region[idxBED])) continue;
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
                if(rv < 0) {
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
                if((plp[0]+i)->qpos+1 > meths[strand-1]->l) meths[strand-1]->l = (plp[0]+i)->qpos+1;
            }
        }
    }

    //Report some output
    if(SVG) makeSVGs(opref, meths, config->keepCpG + 2*config->keepCHG + 4*config->keepCHH);
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

void mbias_usage() {
    fprintf(stderr, "\nUsage: PileOMeth mbias [OPTIONS] <ref.fa> <sorted_alignments.bam> <output.prefix>\n");
    fprintf(stderr,
"\n"
"Options:\n"
" -q INT           Minimum MAPQ threshold to include an alignment (default 10)\n"
" -p INT           Minimum Phred threshold to include a base (default 5). This\n"
"                  must be >0.\n"
" -D INT           Maximum per-base depth (default 2000)\n"
" -r STR           Region string in which to extract methylation\n"
" -l FILE          A BED file listing regions for inclusion. Note that unlike\n"
"                  samtools mpileup, this option will utilize the strand column\n"
"                  (column 6) if present. Thus, if a region has a '+' in this\n"
"                  column, then only metrics from the top strand will be\n"
"                  output. Note that the -r option can be used to limit the\n"
"                  regions of -l.\n"
" --keepDupes      By default, any alignment marked as a duplicate is ignored.\n"
"                  This option causes them to be incorporated.\n"
" --keepSingleton  By default, if only one read in a pair aligns (a singleton)\n"
"                  then it's ignored.\n"
" --keepDiscordant By default, paired-end alignments with the properly-paired bit\n"
"                  unset in the FLAG field are ignored. Note that the definition\n"
"                  of concordant and discordant is based on your aligner\n"
"                  settings.\n"
" -F, --ignoreFlags    By deault, any alignment marked as secondary (bit 0x100),\n"
"                  failing QC (bit 0x200), a PCR/optical duplicate (0x400) or\n"
"                  supplemental (0x800) is ignored. This equates to a value of\n"
"                  0xF00 or 3840 in decimal. If you would like to change that,\n"
"                  you can specify a new value here.\n"
"                  ignored. Specifying this causes them to be included.\n"
" -R, --requireFlags   Require each alignment to have all bits in this value\n"
"                  present, or else the alignment is ignored. This is equivalent\n"
"                  to the -f option in samtools. The default is 0, which\n"
"                  includes all alignments.\n"
" --txt            Output tab separated metrics to the screen. These can be\n"
"                  imported into R or another program for manual plotting and\n"
"                  analysis.\n"
" --noSVG          Don't produce the SVG files. This option implies --txt. Note\n"
"                  that an output prefix is no longer required with this option.\n"
" --noCpG          Do not output CpG methylation metrics\n"
" --CHG            Output CHG methylation metrics\n"
" --CHH            Output CHH methylation metrics\n"
" --version        Print version and the quit\n");
}

int mbias_main(int argc, char *argv[]) {
    char *opref = NULL;
    int c, i, SVG = 1, txt = 0;
    Config config;

    //Defaults
    config.keepCpG = 1; config.keepCHG = 0; config.keepCHH = 0;
    config.minMapq = 10; config.minPhred = 5; config.keepDupes = 0;
    config.keepSingleton = 0, config.keepDiscordant = 0;
    config.maxDepth = 2000;
    config.fai = NULL;
    config.fp = NULL;
    config.bai = NULL;
    config.reg = NULL;
    config.bedName = NULL;
    config.bed = NULL;
    config.ignoreFlags = 0xF00;
    config.requireFlags = 0;
    for(i=0; i<16; i++) config.bounds[i] = 0;

    static struct option lopts[] = {
        {"noCpG",        0, NULL,   1},
        {"CHG",          0, NULL,   2},
        {"CHH",          0, NULL,   3},
        {"keepDupes",    0, NULL,   4},
        {"keepSingleton",  0, NULL, 5},
        {"keepDiscordant", 0, NULL, 6},
        {"txt",          0, NULL,   7},
        {"noSVG",        0, NULL,   8},
        {"ignoreFlags",  1, NULL, 'F'},
        {"requireFlags", 1, NULL, 'R'},
        {"help",         0, NULL, 'h'},
        {"version",      0, NULL, 'v'},
        {0,              0, NULL,   0}
    };
    while((c = getopt_long(argc, argv, "hvq:p:r:l:D:F:", lopts,NULL)) >= 0) {
        switch(c) {
        case 'h' :
            mbias_usage();
            return 0;
        case 'v' :
            print_version();
            return 0;
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
        case 7 :
            txt = 1;
            break;
        case 8 :
            SVG = 0;
            txt = 1;
            break;
        case 'F' :
            config.ignoreFlags = atoi(optarg);
            break;
        case 'R' :
            config.requireFlags = atoi(optarg);
            break;
        case 'q' :
            config.minMapq = atoi(optarg);
            break;
        case 'p' :
            config.minPhred = atoi(optarg);
            break;
        default :
            fprintf(stderr, "Invalid option '%c'\n", c);
            mbias_usage();
            return 1;
        }
    }

    if(argc == 1) {
        mbias_usage();
        return 0;
    }
    if((SVG && argc-optind != 3) || (!SVG && argc-optind < 2)) {
        fprintf(stderr, "You must supply a reference genome in fasta format, an input BAM file, and an output prefix!!!\n");
        mbias_usage();
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
        mbias_usage();
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
    extractMBias(&config, opref, SVG, txt);

    //Close things up
    hts_close(config.fp);
    fai_destroy(config.fai);
    hts_idx_destroy(config.bai);
    if(config.reg) free(config.reg);

    return 0;
}
