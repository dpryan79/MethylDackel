#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>
#include <pthread.h>
#include "MethylDackel.h"

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

strandMeth *mergeStrandMeth(strandMeth *target, strandMeth *source) {
    int32_t i;
    if(source->l == 0) return target;
    if(target->m < source->l) target = growStrandMeth(target, source->m);
    if(target->l < source->l) target->l = source->l;

    for(i=0; i<source->l; i++) {
        target->unmeth1[i] += source->unmeth1[i];
        target->meth1[i] += source->meth1[i];
        target->unmeth2[i] += source->unmeth2[i];
        target->meth2[i] += source->meth2[i];
    }
    return target;
}

void *extractMBias(void *foo) {
    Config *config = (Config*) foo;
    bam_hdr_t *hdr;
    bam_mplp_t iter;
    int ret, tid, pos = 0, i, seqlen, rv, o = 0;
    int32_t bedIdx = 0;
    int strand;
    int n_plp; //This will need to be modified for multiple input files
    const bam_pileup1_t **plp = NULL;
    char *seq = NULL, base;
    mplp_data *data = NULL;
    strandMeth **meths = malloc(4*sizeof(strandMeth*));
    assert(meths);
    uint32_t localPos = 0, localEnd = 0, localTid = 0;
    faidx_t *fai;
    hts_idx_t *bai;
    htsFile *fp;

    for(i=0; i<4; i++) {
        meths[i] = calloc(1, sizeof(strandMeth));
        assert(meths[i]);
    }

    //Open the files
    if((fai = fai_load(config->FastaName)) == NULL) {
        fprintf(stderr, "Couldn't open the index for %s!\n", config->FastaName);
        return NULL;
    }
    if((fp = hts_open(config->BAMName, "rb")) == NULL) {
        fprintf(stderr, "Couldn't open %s for reading!\n", config->BAMName);
        return NULL;
    }
    if((bai = sam_index_load(fp, config->BAMName)) == NULL) {
        fprintf(stderr, "Couldn't load the index for %s\n", config->BAMName);
        return NULL;
    }
    hdr = sam_hdr_read(fp);

    data = calloc(1,sizeof(mplp_data));
    if(data == NULL) {
        fprintf(stderr, "Couldn't allocate space for the data structure in extractCalls()!\n");
        return NULL;
    }
    data->config = config;
    data->hdr = hdr;
    data->fp = fp;
    data->bedIdx = bedIdx;

    plp = calloc(1, sizeof(bam_pileup1_t *)); //This will have to be modified for multiple input files
    if(plp == NULL) {
        fprintf(stderr, "Couldn't allocate space for the plp structure in extractCalls()!\n");
        return NULL;
    }

    while(1) {
        //Lock and unlock the mutex so we can get/update the tid/position
        pthread_mutex_lock(&positionMutex);
        localTid = globalTid;
        localPos = globalPos;
        localEnd = localPos + config->chunkSize;
        if(localTid >= hdr->n_targets) {
            pthread_mutex_unlock(&positionMutex);
            break;
        }
        if(globalEnd && localEnd > globalEnd) localEnd = globalEnd;
        adjustBounds(config, hdr, fai, &localTid, &localPos, &localEnd);
        globalPos = localEnd;
        if(globalEnd > 0 && globalPos >= globalEnd) {
            //If we've specified a region, then break once we're outside of it
            globalTid = (uint32_t) -1;
        }
        if(localTid < hdr->n_targets && globalTid != (uint32_t) -1) {
            if(globalPos >= hdr->target_len[localTid]) {
                localEnd = hdr->target_len[localTid];
                globalTid++;
                globalPos = 0;
            }
        }
        pthread_mutex_unlock(&positionMutex);

        //If we have a BED file, then jump to the first overlapping region
        if(config->bed) {
            if(spanOverlapsBED(localTid, localPos, localEnd, config->bed, &bedIdx) != 1) continue;
        }

        //Break out of the loop if finished
        if(localTid >= hdr->n_targets) break; //Finish looping
        if(globalEnd && localPos >= globalEnd) break;
        data->iter = sam_itr_queryi(bai, localTid, localPos, localEnd);

        seq = faidx_fetch_seq(fai, hdr->target_name[localTid], localPos, localEnd, &seqlen);
        if(seqlen < 0) {
            fprintf(stderr, "faidx_fetch_seq returned %i while trying to fetch the sequence for tid %s:%"PRIu32"-%"PRIu32"!\n",\
                seqlen, hdr->target_name[localTid], localPos, localEnd);
            fprintf(stderr, "Note that the output will be truncated!\n");
            return NULL;
        }
        data->seq = seq;
        data->lseq = seqlen;
        data->offset = localPos;

        //Start the pileup
        iter = bam_mplp_init(1, filter_func, (void **) &data);
//        bam_mplp_init_overlaps(iter); //This is included in extract but excluded here. The main benefit to exclusion is that you can more accurately gauge overlapping regions.
        bam_mplp_set_maxcnt(iter, INT_MAX);
        while((ret = bam_mplp_auto(iter, &tid, &pos, &n_plp, plp)) > 0) {
            if(pos < localPos || pos >= localEnd) continue; // out of the region requested

            if(config->bed) { //Handle -l
                while((o = posOverlapsBED(tid, pos, config->bed, bedIdx)) == -1) bedIdx++;
                if(o == 0) continue; //Wrong strand
            }

            if(isCpG(seq, pos-localPos, seqlen)) {
                if(!config->keepCpG) continue;
            } else if(isCHG(seq, pos-localPos, seqlen)) {
                if(!config->keepCHG) continue;
            } else if(isCHH(seq, pos-localPos, seqlen)) {
                if(!config->keepCHH) continue;
            } else {
                continue;
            }

            base = *(seq+pos-localPos);
            for(i=0; i<n_plp; i++) {
                if(plp[0][i].is_del) continue;
                if(plp[0][i].is_refskip) continue;
                if(config->bed) if(!readStrandOverlapsBED(plp[0][i].b, config->bed->region[bedIdx])) continue;
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
        hts_itr_destroy(data->iter);
        free(seq);
        bam_mplp_destroy(iter);
    }

    //Clean up
    bam_hdr_destroy(hdr);
    fai_destroy(fai);
    hts_close(fp);
    hts_idx_destroy(bai);
    free(data);
    free(plp);

    return meths;
}

void mbias_usage() {
    fprintf(stderr, "\nUsage: MethylDackel mbias [OPTIONS] <ref.fa> <sorted_alignments.bam> <output.prefix>\n");
    fprintf(stderr,
"\n"
"Options:\n"
" -q INT           Minimum MAPQ threshold to include an alignment (default 10)\n"
" -p INT           Minimum Phred threshold to include a base (default 5). This\n"
"                  must be >0.\n"
" -D INT           Maximum per-base depth (default 2000)\n"
" --minIsize INT   Minimum insert size / fragment length. May be useful for\n"
"                  applications such as cell-free DNA. Inclusive.\n"
" --maxIsize INT   Maximum insert size / fragment length. Also inclusive.\n"
" -r STR           Region string in which to extract methylation\n"
" -l FILE          A BED file listing regions for inclusion.\n"
" --keepStrand     If a BED file is specified, then this option will cause the\n"
"                  strand column (column 6) to be utilized, if present. Thus, if\n"
"                  a region has a '+' in this column, then only metrics from the\n"
"                  top strand will be output. Note that the -r option can be used\n"
"                  to limit the regions of -l.\n"
" -@ nThreads      The number of threads to use, the default 1\n"
" --chunkSize INT  The size of the genome processed by a single thread at a time.\n"
"                  The default is 1000000 bases. This value MUST be at least 1.\n"
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
" --ignoreNH       Ignore NH auxiliary tags. By default, if an NH tag is present\n"
"                  and its value is >1 then an entry is ignored as a\n"
"                  multimapper.\n"
" --minConversionEfficiency  The minimum non-CpG conversion efficiency observed\n"
"                  in a read to include it in the output. The default is 0.0 and\n"
"                  the maximum is 1.0 (100%% conversion). You are strongly\n"
"                  encouraged to NOT use this option without an EXTREMELY\n"
"                  compelling reason!\n"
" --txt            Output tab separated metrics to the screen. These can be\n"
"                  imported into R or another program for manual plotting and\n"
"                  analysis. Note that coordinates are 1-based.\n"
" --noSVG          Don't produce the SVG files. This option implies --txt. Note\n"
"                  that an output prefix is no longer required with this option.\n"
" --noCpG          Do not output CpG methylation metrics\n"
" --CHG            Output CHG methylation metrics\n"
" --CHH            Output CHH methylation metrics\n"
" --nOT INT,INT,INT,INT Inclusion bound for methylation calls from reads/pairs\n"
"                  originating from the original top strand. Each integer\n"
"                  represents a 1-based position from the end of a read. For\n"
"                  example \"--nOT A,B,C,D\" translates to, \"Include calls from\n"
"                  position A through the Bth read from the end on read #1 and\n"
"                  Cth through the Dth from the end base on read #2\". In other\n"
"                  words \"--nOT 5,10,0,0\" for a 100 base long read would result\n"
"                  in bases 5 through 90 being used. If a 0 is used in any\n"
"                  position then that is translated to mean start/end of the\n"
"                  alignment, as appropriate. For example, --nOT 5,0,0,0 would\n"
"                  include all but the first 4 bases on read #1.\n"
" --nOB INT,INT,INT,INT\n"
" --nCTOT INT,INT,INT,INT\n"
" --nCTOB INT,INT,INT,INT As with --nOT, but for the original bottom, complementary\n"
"                  to the original top, and complementary to the original bottom\n"
"                  strands, respectively.\n"
" --fivePrime INT  Alternative trimming option to --OT / --nOT. Trimming based on\n"
"                  fragment ends rather than read ends. Experimental, and is not\n"
"                  accurate in cases where trim length is greater than read length\n"
" --threePrime INT\n"
" --version        Print version and the quit\n");
}

int mbias_main(int argc, char *argv[]) {
    char *opref = NULL;
    int c, i, j, SVG = 1, txt = 0, keepStrand = 0;
    strandMeth *meths[4], **threadout = NULL;
    Config config;
    bam_hdr_t *hdr = NULL;

    //Defaults
    config.keepCpG = 1; config.keepCHG = 0; config.keepCHH = 0;
    config.minMapq = 10; config.minPhred = 5; config.keepDupes = 0;
    config.keepSingleton = 0, config.keepDiscordant = 0;
    config.filterMappability = 0, config.ignoreNH = 0;
    
    config.fp = NULL;
    config.bai = NULL;
    config.reg = NULL;
    config.bedName = NULL;
    config.bed = NULL;
    config.ignoreFlags = 0xF00;
    config.requireFlags = 0;
    config.nThreads = 1;
    config.chunkSize = 1000000;
    config.minConversionEfficiency = 0.0;
    for(i=0; i<16; i++) config.bounds[i] = 0;
    for(i=0; i<16; i++) config.absoluteBounds[i] = 0;
    config.fivePrime = 0;
    config.threePrime = 0;
    config.minIsize = 0;
    config.maxIsize = 0;

    static struct option lopts[] = {
        {"noCpG",        0, NULL,   1},
        {"CHG",          0, NULL,   2},
        {"CHH",          0, NULL,   3},
        {"keepDupes",    0, NULL,   4},
        {"keepSingleton",  0, NULL, 5},
        {"keepDiscordant", 0, NULL, 6},
        {"txt",          0, NULL,   7},
        {"noSVG",        0, NULL,   8},
        {"nOT",          1, NULL,   9},
        {"nOB",          1, NULL,  10},
        {"nCTOT",        1, NULL,  11},
        {"nCTOB",        1, NULL,  12},
        {"chunkSize",    1, NULL,  13},
        {"keepStrand",   0, NULL,  14},
        {"minConversionEfficiency", 1, NULL, 15},
        {"ignoreNH",     0, NULL,  16},
        {"fivePrime",  1, NULL, 17},
        {"threePrime", 1, NULL, 18},
        {"minIsize",  1, NULL, 19},
        {"maxIsize", 1, NULL, 20},
        {"ignoreFlags",  1, NULL, 'F'},
        {"requireFlags", 1, NULL, 'R'},
        {"help",         0, NULL, 'h'},
        {"version",      0, NULL, 'v'},
        {0,              0, NULL,   0}
    };
    while((c = getopt_long(argc, argv, "hvq:p:r:l:D:F:@:", lopts,NULL)) >= 0) {
        switch(c) {
        case 'h' :
            mbias_usage();
            return 0;
        case 'v' :
            print_version();
            return 0;
        case 'D' :
            // This is now set to INT_MAX, it was --maxDepth
            break;
        case 'r':
            config.reg = optarg;
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
        case 9 :
            parseBounds(optarg, config.absoluteBounds, 0);
            break;
        case 10 :
            parseBounds(optarg, config.absoluteBounds, 1);
            break;
        case 11 :
            parseBounds(optarg, config.absoluteBounds, 2);
            break;
        case 12 :
            parseBounds(optarg, config.absoluteBounds, 3);
            break;
        case 13:
            config.chunkSize = strtoul(optarg, NULL, 10);
            if(config.chunkSize < 1) {
                fprintf(stderr, "Error: The chunk size must be at least 1!\n");
                return 1;
            }
            break;
        case 14:
            keepStrand = 1;
            break;
        case 15:
            config.minConversionEfficiency = atof(optarg);
            break;
        case 16:
            config.ignoreNH = 1;
            break;
        case 17:
            config.fivePrime = atoi(optarg);
            break;
        case 18:
            config.threePrime = atoi(optarg);
            break;
        case 19:
            config.minIsize = atoi(optarg);
            break;
        case 20:
            config.maxIsize = atoi(optarg);
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
        case '@':
            config.nThreads = atoi(optarg);
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
    if(config.fivePrime < 0) {
        fprintf(stderr, "--fivePrime %i is invalid. Resetting to 0, which is the lowest possible value.\n", config.fivePrime);
        config.fivePrime = 0;
    }
    if(config.threePrime < 0) {
        fprintf(stderr, "--threePrime %i is invalid. Resetting to 0, which is the lowest possible value.\n", config.threePrime);
        config.threePrime = 0;
    }
    if(config.minIsize < 0) {
        fprintf(stderr, "--minIsize %i is invalid. Resetting to 0, which is the default value.\n", config.minIsize);
        config.minIsize = 0;
    }
    if(config.maxIsize < 0) {
        fprintf(stderr, "--maxIsize %i is invalid. Resetting to 0, which is the default value.\n", config.maxIsize);
        config.maxIsize = 0;
    }

    //Is there still a metric to output?
    if(!(config.keepCpG + config.keepCHG + config.keepCHH)) {
        fprintf(stderr, "You haven't specified any metrics to output!\nEither don't use the --noCpG option or specify --CHG and/or --CHH.\n");
        return -1;
    }

    //Open the files
    config.FastaName = argv[optind];
    config.BAMName = argv[optind+1];
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

    //parse the region, if needed
    if(config.reg) {
        const char *foo;
        char *bar = NULL;
        int s, e;
        hdr = sam_hdr_read(config.fp);
        foo = hts_parse_reg(config.reg, &s, &e);
        if(foo == NULL) {
            fprintf(stderr, "Could not parse the specified region!\n");
            return -4;
        }
        bar = malloc(foo - config.reg + 1);
        if(bar == NULL) {
            fprintf(stderr, "Could not allocate temporary space for parsing the requested region!\n");
            return -5;
        }
        strncpy(bar, config.reg, foo - config.reg);
        bar[foo - config.reg] = 0;
        globalTid = bam_name2id(hdr, bar);
        if(globalTid == (uint32_t) -1) {
            fprintf(stderr, "%s did not match a known chromosome/contig name!\n", config.reg);
            return -6;
        }
        if(s>0) globalPos = s;
        if(e>0) globalEnd = e;
        if(globalEnd > hdr->target_len[globalTid]) globalEnd = hdr->target_len[globalTid];
        free(bar);
        if(!config.bedName) bam_hdr_destroy(hdr);
    }
    if(config.bedName) {
        if(!hdr) hdr = sam_hdr_read(config.fp);
        config.bed = parseBED(config.bedName, hdr, keepStrand);
        bam_hdr_destroy(hdr);
        if(!config.bed) {
            fprintf(stderr, "There was an error while reading in your BED file!\n");
            return 1;
        }
    }

    for(i=0; i<4; i++) {
        meths[i] = calloc(1, sizeof(strandMeth));
        assert(meths[i]);
    }

    //Run the pileup
    pthread_mutex_init(&positionMutex, NULL);
    pthread_t *threads = calloc(config.nThreads, sizeof(pthread_t));
    for(i=0; i < config.nThreads; i++) pthread_create(threads+i, NULL, &extractMBias, &config);
    for(i=0; i < config.nThreads; i++) {
        pthread_join(threads[i], (void**) &threadout);
        for(j=0; j<4; j++) {
            meths[j] = mergeStrandMeth(meths[j], threadout[j]);
            if(threadout[j]->meth1) free(threadout[j]->meth1);
            if(threadout[j]->unmeth1) free(threadout[j]->unmeth1);
            if(threadout[j]->meth2) free(threadout[j]->meth2);
            if(threadout[j]->unmeth2) free(threadout[j]->unmeth2);
            free(threadout[j]);
        }
        free(threadout);
    }

    free(threads);
    pthread_mutex_destroy(&positionMutex);

    //Report some output
    if(SVG) makeSVGs(opref, meths, config.keepCpG + 2*config.keepCHG + 4*config.keepCHH);
    if(txt) makeTXT(meths);

    //Close things up
    hts_close(config.fp);
    hts_idx_destroy(config.bai);
    for(i=0; i<4; i++) {
        if(meths[i]->meth1) free(meths[i]->meth1);
        if(meths[i]->unmeth1) free(meths[i]->unmeth1);
        if(meths[i]->meth2) free(meths[i]->meth2);
        if(meths[i]->unmeth2) free(meths[i]->unmeth2);
        free(meths[i]);
    }

    return 0;
}
