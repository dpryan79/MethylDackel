#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>
#include <limits.h>
#include <pthread.h>
#include "MethylDackel.h"

void print_version(void);

void addRead(kstring_t *os, bam1_t *b, bam_hdr_t *hdr, uint32_t nmethyl, uint32_t nunmethyl) {
    char str[10000]; // I don't really like hardcoding it, but given the probability that it ever won't suffice...

    if(nmethyl + nunmethyl > 0) {
        snprintf(str, 10000, "%s\t%s\t%i\t%f\t%"PRIu32"\n",
            bam_get_qname(b),
            hdr->target_name[b->core.tid],
            b->core.pos,
            100. * ((double) nmethyl)/(nmethyl+nunmethyl),
            nmethyl + nunmethyl);
    } else {
        snprintf(str, 10000, "%s\t%s\t%i\t0.0\t%"PRIu32"\n",
            bam_get_qname(b),
            hdr->target_name[b->core.tid],
            b->core.pos,
            nmethyl + nunmethyl);
    }

    kputs(str, os);
}

void processRead(Config *config, bam1_t *b, char *seq, uint32_t sequenceStart, int seqLen, uint32_t *nmethyl, uint32_t *nunmethyl) {
    uint32_t readPosition = 0;
    uint32_t mappedPosition = b->core.pos;
    int cigarOPNumber = 0;
    int cigarOPOffset = 0;
    uint32_t *CIGAR = bam_get_cigar(b);
    uint8_t *readSeq = bam_get_seq(b); //get_get_qual() too?
    int strand = getStrand(b);
    int cigarOPType;
    int direction;
    int base;

    while(readPosition < b->core.l_qseq && cigarOPOffset < b->core.n_cigar) {
        if(cigarOPOffset >= bam_cigar_oplen(CIGAR[cigarOPNumber])) {
            cigarOPOffset = 0;
            cigarOPNumber++;
        }
        cigarOPType = bam_cigar_type(CIGAR[cigarOPNumber]);
        if(cigarOPType & 2) { //not ISHPB
            if(cigarOPType & 1) { //M=X
                direction = isCpG(seq, mappedPosition - sequenceStart, seqLen);
                if(direction) {
                    base = bam_seqi(readSeq, readPosition);  // Filtering by quality goes here
                    if(direction == 1 && (strand & 1) == 1) { // C & OT/CTOT
                        if(base == 2) (*nmethyl)++;  //C
                        else if(base == 8) (*nunmethyl)++; //T
                    } else if(direction == -1 && (strand & 1) == 0) { // G & OB/CTOB
                        if(base == 4) (*nmethyl)++;  //G
                        else if(base == 1) (*nunmethyl)++; //A
                    }
                }
                mappedPosition++;
                readPosition++;
                cigarOPOffset++;
            } else { //DN
                mappedPosition += bam_cigar_oplen(CIGAR[cigarOPNumber++]);
                cigarOPOffset = 0;
                continue;
            }
        } else if(cigarOPType & 1) { // IS
            readPosition += bam_cigar_oplen(CIGAR[cigarOPNumber++]);
            cigarOPOffset = 0;
            continue;
        } else { // HPB Note that B is not handled properly, but it doesn't currently exist in the wild
            cigarOPOffset = 0;
            cigarOPNumber++;
            continue;
        }
    }
}

void *perReadMetrics(void *foo) {
    Config *config = (Config*) foo;
    bam_hdr_t *hdr;
    int32_t bedIdx = 0;
    int seqlen;
    uint32_t nmethyl = 0, nunmethyl = 0;
    uint32_t localBin = 0, localPos = 0, localEnd = 0, localTid = 0, localPos2 = 0;
    char *seq = NULL;
    kstring_t *os = NULL;
    faidx_t *fai;
    hts_idx_t *bai;
    htsFile *fp;
    bam1_t *b = bam_init1();
    hts_itr_t *iter;

    os = calloc(1, sizeof(kstring_t));
    if(!os) {
        fprintf(stderr, "Couldn't allocate space the for kstring_t structure in perRead!\n");
        return NULL;
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

    while(1) {
        //Lock and unlock the mutex so we can get/update the tid/position
        pthread_mutex_lock(&positionMutex);
        localBin = bin++;
        localTid = globalTid;
        localPos = globalPos;
        localEnd = localPos + config->chunkSize;
        if(localTid >= hdr->n_targets) {
            pthread_mutex_unlock(&positionMutex);
            break;
        }
        if(globalEnd && localEnd > globalEnd) localEnd = globalEnd;
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

        //If we have a BED file, then jump to the first overlapping region:
        if(config->bed) {
            if(spanOverlapsBED(localTid, localPos, localEnd, config->bed, &bedIdx) != 1) {
                //Set the bin as written and loop
                while(1) {
                    pthread_mutex_lock(&outputMutex);
                    if(outputBin != localBin) {
                        pthread_mutex_unlock(&outputMutex);
                        continue;
                    }
                    outputBin++;
                    break;
                }
                pthread_mutex_unlock(&outputMutex);
                continue;
            }
        }
        localPos2 = 0;
        if(localPos > 1) {
            localPos2 = localPos - 2;
        }

        //Break out of the loop if finished
        if(localTid >= hdr->n_targets) break; //Finish looping
        if(globalEnd && localPos >= globalEnd) break;
        iter = sam_itr_queryi(bai, localTid, localPos, localEnd);

        //Fetch the sequence and then iterate over the reads. There's a 10kb buffer on the end
        seq = faidx_fetch_seq(fai, hdr->target_name[localTid], localPos2, localEnd+10000, &seqlen);
        while(sam_itr_next(fp, iter, b) >= 0) {
            if(b->core.pos < localPos2) continue;
            nmethyl = 0, nunmethyl = 0;
            processRead(config, b, seq, localPos2, seqlen, &nmethyl, &nunmethyl);
            addRead(os, b, hdr, nmethyl, nunmethyl);
        }
        sam_itr_destroy(iter);
        free(seq);

        //Write output
        while(1) {
            pthread_mutex_lock(&outputMutex);
            if(outputBin != localBin) {
                pthread_mutex_unlock(&outputMutex);
                continue;
            }
            fputs(os->s, config->output_fp[0]);
            os->l = 0;
            outputBin++;
            pthread_mutex_unlock(&outputMutex);
            break;
        }
    }

    if(os->s) free(os->s);
    free(os);
    bam_hdr_destroy(hdr);
    bam_destroy1(b);
    fai_destroy(fai);
    hts_close(fp);
    hts_idx_destroy(bai);
    return NULL;
}

void perRead_usage() {
    fprintf(stderr, "\nUsage: MethylDackel perRead [OPTIONS] <ref.fa> <input>\n");
    fprintf(stderr,
"\n"
"This program will compute the average CpG methylation level of a given read.\n"
"The output is a tab-separated file with the following columns:\n"
" - read name\n"
" - chromosome\n"
" - position\n"
" - CpG methylation (%%)\n"
" - number of informative bases\n"
"\n"
"Arguments:\n"
"  ref.fa    Reference genome in fasta format. This must be indexed with\n"
"            samtools faidx\n"
"  input     An input BAM or CRAM file. This MUST be sorted and should be indexed.\n"
"\nOptions:\n"
"  -o STR    Output file name [stdout]\n"
" -@ nThreads      The number of threads to use, the default 1\n"
" --chunkSize INT  The size of the genome processed by a single thread at a time.\n"
"                  The default is 1000000 bases. This value MUST be at least 1.\n"
"  --version Printer version and quit\n"
"\n"
"Note that this program will produce incorrect values for alignments spanning\n"
"more than 10kb.\n");
}

int perRead_main(int argc, char *argv[]) {
    Config config;
    FILE *ofile = stdout;
    int i, keepStrand = 0;
    faidx_t *fai;
    bam_hdr_t *hdr = NULL;
    char c;

    //Defaults
    config.keepCpG = 1; config.keepCHG = 0; config.keepCHH = 0;
    config.minMapq = 0; config.minPhred = 0; config.keepDupes = 0;
    config.keepSingleton = 0, config.keepDiscordant = 0;
    config.fp = NULL;
    config.bai = NULL;
    config.reg = NULL;
    config.bedName = NULL;
    config.bed = NULL;
    config.ignoreFlags = 0;
    config.requireFlags = 0;
    config.nThreads = 1;
    config.chunkSize = 1000000;

    static struct option lopts[] = {
        {"help",    0, NULL, 'h'},
        {"version", 0, NULL, 'v'},
        {"chunkSize",    1, NULL,  19},
        {0,         0, NULL,   0}
    };
    //Add filtering options
    //BED file support
    //region support
    //stdout vs. file name
    while((c = getopt_long(argc, argv, "hvo:", lopts, NULL)) >= 0) {
        switch(c) {
        case 'h' :
            perRead_usage();
            return 0;
        case 'v' :
            print_version();
            return 0;
        case 'o' :
            if((ofile = fopen(optarg, "w")) == NULL) {
                fprintf(stderr, "Couldn't open %s for writing\n", optarg);
                return 2;
            }
            break;
        case '@':
            config.nThreads = atoi(optarg);
            break;
        case 19:
            config.chunkSize = strtoul(optarg, NULL, 10);
            if(config.chunkSize < 1) {
                fprintf(stderr, "Error: The chunk size must be at least 1!\n");
                return 1;
            }
            break;
        default :
            fprintf(stderr, "Invalid option '%c'\n", c);
            perRead_usage();
            return 1;
        }
    }

    if(argc == 1) {
        perRead_usage();
        return 0;
    }
    if(argc - optind != 2) {
        fprintf(stderr, "You must supply a reference genome in fasta format and a BAM or CRAM file\n");
        perRead_usage();
        return -1;
    }
    if((fai = fai_load(argv[optind])) == NULL) {
        fprintf(stderr, "Couldn't open the index for %s!\n", argv[optind]);
        perRead_usage();
        return -2;
    }
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

    //Output file
    config.output_fp = malloc(sizeof(FILE *));
    config.output_fp[0] = ofile;

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

    //Spin up the threads
    pthread_mutex_init(&positionMutex, NULL);
    pthread_t *threads = calloc(config.nThreads, sizeof(pthread_t));
    for(i=0; i < config.nThreads; i++) pthread_create(threads+i, NULL, &perReadMetrics, &config);
    for(i=0; i < config.nThreads; i++) pthread_join(threads[i], NULL);
    free(threads);
    pthread_mutex_destroy(&positionMutex);

    //Close things up
    hts_close(config.fp);
    if(ofile != stdout) fclose(ofile);
    free(config.output_fp);
    fai_destroy(fai);
    hts_idx_destroy(config.bai);
    if(config.bed) destroyBED(config.bed);

    return 0;
}
