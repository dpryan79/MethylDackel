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

    // parsing mate cigar string
    char *c_string = bam_aux2Z(bam_aux_get(b, "MC"));
    char *end;
    uint32_t *mate_cigar = NULL;
    size_t m = 0;
    int n_cigar_mate;
    n_cigar_mate = sam_parse_cigar(c_string, &end, &mate_cigar, &m);
    int uncovered = abs(b->core.isize) - bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b)) - bam_cigar2rlen(n_cigar_mate, mate_cigar);
    
    snprintf(str, 10000, "%s\t%s\t%"PRId64"\t%"PRIu32"\t%"PRIu32"\t%"PRId32"\t%"PRId64"\t%"PRIu16"\t%"PRId64"\t%"PRId64"\t%s\n",
        bam_get_qname(b),
        hdr->target_name[b->core.tid],
        b->core.pos,
        nmethyl,
        nunmethyl,
        (uncovered > 0) ? uncovered : 0,
        b->core.isize,
        b->core.flag,
        bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b)),
        bam_cigar2rlen(n_cigar_mate, mate_cigar),
        c_string
    );
    free(mate_cigar);
    kputs(str, os);
}
/*
int calcFragmentBound(int readPosWanted, int n_cigar, uint32_t *CIGAR, int mapPos) {
    int readPos = 0;
    int cigarPos = 0;
    int cigarOPType;
    int qlen = bam_cigar2qlen(n_cigar, CIGAR);
    
    if(readPosWanted > 0) {
        while (readPos < qlen && cigarPos < n_cigar) {
            cigarOPType = bam_cigar_type(CIGAR[cigarPos]);
            if(cigarOPType & 2) {
                if(cigarOPType & 1) { // consumes both query and ref
                    if((readPos + bam_cigar_oplen(CIGAR[cigarPos])) >= readPosWanted) {
                        mapPos = mapPos + readPosWanted - readPos; // move mapPos to the location where readPosWanted is
                        break;
                    } else {
                        readPos += bam_cigar_oplen(CIGAR[cigarPos]);
                        mapPos += bam_cigar_oplen(CIGAR[cigarPos++]);
                        continue;
                    }
                } else { // consumes only ref but not query
                    mapPos += bam_cigar_oplen(CIGAR[cigarPos++]);
                    continue;
                }
            } else if(cigarOPType & 1) { // consumes only query but not ref
                readPos += bam_cigar_oplen(CIGAR[cigarPos++]);
                if(readPos >= readPosWanted) {
                    break; // mapPos does not have to be moved
                } else {
                    continue;
                }
            } else { // consumes neither query or ref
                cigarPos++;
                continue;
            }
        }
    }
    return(mapPos);
}
*/
void processRead(Config *config, bam1_t *b, char *seq, uint32_t sequenceStart, int seqLen, uint32_t *nmethyl, uint32_t *nunmethyl) {
    uint32_t readPosition = 0;
    uint32_t mappedPosition = b->core.pos;
    int cigarOPNumber = 0;
    int cigarOPOffset = 0;
    uint32_t *CIGAR = bam_get_cigar(b);
    uint8_t *readSeq = bam_get_seq(b);
    uint8_t *readQual = bam_get_qual(b);
    int strand = getStrand(b);
    int cigarOPType;
    int direction;
    int base;

    // parsing mate cigar string
    char *c_string = bam_aux2Z(bam_aux_get(b, "MC"));
    char *end;
    uint32_t *mate_cigar = NULL;
    size_t m = 0;
    int n_cigar_mate;
    n_cigar_mate = sam_parse_cigar(c_string, &end, &mate_cigar, &m);
    uint32_t matePosition = b->core.mpos + bam_cigar2rlen(n_cigar_mate, mate_cigar);
/*    int mateQlen = bam_cigar2qlen(n_cigar_mate, mate_cigar);

    // 5' and 3' trimming, assumes paired end sequencing
    int lb, rb;
    
    // Left and right bounds based on desired 5' and 3' trimming
    // Original version based on ref position, which is more straightforward
    // Can be further simplified but it's hard enough to read already
    if(b->core.flag & 0x3) { // Read properly paired
      if(b->core.flag & 0x40) { // Read 1 (ie contains 5')
          if(b->core.flag & 0x20) { // OT
              lb = calcFragmentBound(config->fivePrime, b->core.n_cigar, bam_get_cigar(b), b->core.pos);
              rb = calcFragmentBound(mateQlen - config->threePrime, n_cigar_mate, mate_cigar, b->core.mpos);
          } else { // OB
              lb = calcFragmentBound(config->threePrime, n_cigar_mate, mate_cigar, b->core.mpos);
              rb = calcFragmentBound(b->core.l_qseq - config->fivePrime, b->core.n_cigar, bam_get_cigar(b), b->core.pos);
          }
      } else { // Read 2 (ie contains 3')
          if(b->core.flag & 0x20) { // OT
              lb = calcFragmentBound(config->threePrime, b->core.n_cigar, bam_get_cigar(b), b->core.pos);
              rb = calcFragmentBound(mateQlen - config->fivePrime, n_cigar_mate, mate_cigar, b->core.mpos);
          } else { // OB
              lb = calcFragmentBound(config->fivePrime, n_cigar_mate, mate_cigar, b->core.mpos);
              rb = calcFragmentBound(b->core.l_qseq - config->threePrime, b->core.n_cigar, bam_get_cigar(b), b->core.pos);
          }
      }
    } else {
        fprintf(stderr, "Read with qname [[ %s ]] not properly paired and is skipped\n", bam_get_qname(b));
        lb = -1;
        rb = -1;
    }
*/
    while(readPosition < b->core.l_qseq && cigarOPNumber < b->core.n_cigar ) { //&& mappedPosition <= rb) {
        if(cigarOPOffset >= bam_cigar_oplen(CIGAR[cigarOPNumber])) {
            cigarOPOffset = 0;
            cigarOPNumber++;
        }
        cigarOPType = bam_cigar_type(CIGAR[cigarOPNumber]);
        if(cigarOPType & 2) { //not ISHPB
            if(cigarOPType & 1) { //M=X
                // If read is on the reverse strand (ie rightward to its mate), skip overlap with mate
                // Ignorant to phread score; there's a better way of doing this but too cumbersome 
                if(bam_is_rev(b) && (mappedPosition <= matePosition)) {
                    mappedPosition++;
                    readPosition++;
                    cigarOPOffset++;
                    continue;
                }
/*
                // Skip the left bound; right bound is skipped using the while loop conditions
                if(mappedPosition < lb) {
                    mappedPosition++;
                    readPosition++;
                    cigarOPOffset++;
                    continue;
                }
*/
                // Skip poor base calls
                if(readQual[readPosition] < config->minPhred) {
                    mappedPosition++;
                    readPosition++;
                    cigarOPOffset++;
                    continue;
                }

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
    free(mate_cigar);
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
            if(b->core.pos < localPos) continue;
            if(b->core.pos >= localEnd) break;
            nmethyl = 0, nunmethyl = 0;
            if(config->requireFlags && (config->requireFlags & b->core.flag) != config->requireFlags) continue;
            if(config->ignoreFlags && (config->ignoreFlags & b->core.flag) != 0) continue;
            if(b->core.qual < config->minMapq) continue;
            if(config->fivePrime || config->threePrime) b = trimFragmentEnds(b, config->fivePrime, config->threePrime);
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
            if(os->l) fputs(os->s, config->output_fp[0]);
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
" - number of unmodified cytosine (ie 5mC for bisulfite)\n"
" - number of modified cytosine\n"
" - alignment overlap\n"
" - insert size\n"
" - sam flag (for strand & read information)\n"
" - self rlen\n"
" - mate rlen\n"
"\n"
"Arguments:\n"
"  ref.fa    Reference genome in fasta format. This must be indexed with\n"
"            samtools faidx\n"
"  input     An input BAM or CRAM file. This MUST be sorted and should be indexed.\n"
"\nOptions:\n"
" -q INT     Minimum MAPQ threshold to include an alignment (default 10)\n"
" -p INT     Minimum Phred threshold to include a base (default 5). This must\n"
"            be >0.\n"
" -r STR     Region string in which to extract methylation\n"
" -l FILE    A BED file listing regions for inclusion.\n"
" --keepStrand  If a BED file is specified, then this option will cause the\n"
"            strand column (column 6) to be utilized, if present. Thus, if\n"
"            a region has a '+' in this column, then only metrics from the\n"
"            top strand will be output. Note that the -r option can be used\n"
"            to limit the regions of -l.\n"
" -o STR     Output file name [stdout]\n"
" -F, --ignoreFlags    By default, all reads are output. If you would like to\n"
"            ignore certain classes of reads then simply give a value for their\n"
"            flags here. Note that an alignment will be logically anded with this\n"
"            flag, so a single bit overlap will lead to exclusion. The default\n"
"            for this is 0.\n"
" -R, --requireFlags   Require each alignment to have all bits in this value\n"
"            present, or else the alignment is ignored. This is equivalent to the\n"
"            -f option in samtools. The default is 0, which includes all\n"
"            alignments.\n"
" --ignoreNH Ignore NH auxiliary tags. By default, if an NH tag is present\n"
"            and its value is >1 then an entry is ignored as a\n"
"            multimapper.\n"
" -@ INT     The number of threads to use, the default 1\n"
" --chunkSize INT  The size of the genome processed by a single thread at a time.\n"
"            The default is 1000000 bases. This value MUST be at least 1.\n"
" --fivePrime INT  Alternative trimming option to --OT / --nOT. Trimming based on\n"
"                  fragment ends rather than read ends. Experimental, and is not\n"
"                  accurate in cases where trim length is greater than read length\n"
" --threePrime INT\n"
" --minIsize INT  Filter by minimum insert size; inclusive of INT\n"
" --maxIsize INT  Filter by maximum insert size; also inclusive\n"
" --version  Print version and quit\n"
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
    config.minMapq = 10; config.minPhred = 5; config.keepDupes = 0;
    config.keepSingleton = 0, config.keepDiscordant = 0;
    config.ignoreNH = 0;
    config.fp = NULL;
    config.bai = NULL;
    config.reg = NULL;
    config.bedName = NULL;
    config.bed = NULL;
    config.ignoreFlags = 0;
    config.requireFlags = 0;
    config.nThreads = 1;
    config.chunkSize = 1000000;
    config.fivePrime = 0;
    config.threePrime = 0;

    static struct option lopts[] = {
        {"help",    0, NULL, 'h'},
        {"version", 0, NULL, 'v'},
        {"chunkSize",    1, NULL,  19},
        {"keepStrand",   0, NULL,  20},
        {"ignoreFlags",  1, NULL, 'F'},
        {"requireFlags", 1, NULL, 'R'},
        {"fivePrime",  1, NULL, 22},
        {"threePrime", 1, NULL, 23},
        {"minIsize", 1, NULL, 24},
        {"maxIsize", 1, NULL, 25},
        {0,         0, NULL,   0}
    };
    while((c = getopt_long(argc, argv, "hvq:p:o:@:r:l:F:R:", lopts, NULL)) >= 0) {
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
        case 'q':
            config.minMapq = atoi(optarg);
            break;
        case 'p':
            config.minPhred = atoi(optarg);
            break;
        case '@':
            config.nThreads = atoi(optarg);
            break;
        case 'r':
            config.reg = optarg;
            break;
        case 'l':
            config.bedName = optarg;
            break;
        case 'F':
            config.ignoreFlags = atoi(optarg);
            break;
        case 'R':
            config.requireFlags = atoi(optarg);
            break;
        case 19:
            config.chunkSize = strtoul(optarg, NULL, 10);
            if(config.chunkSize < 1) {
                fprintf(stderr, "Error: The chunk size must be at least 1!\n");
                return 1;
            }
            break;
        case 20:
            keepStrand = 1;
            break;
        case 21:
            config.ignoreNH = 1;
            break;
        case 22:
            config.fivePrime = atoi(optarg);
            break;
        case 23:
            config.threePrime = atoi(optarg);
            break;
        case 24:
            config.minIsize = atoi(optarg);
            break;
        case 25:
            config.maxIsize = atoi(optarg);
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
        fprintf(stderr, "--minIsize %i is invalid. Resetting to 0, which will not filter based on min insert size.\n", config.minIsize);
        config.minIsize = 0;
    }
    if(config.maxIsize < 0) {
        fprintf(stderr, "--maxIsize %i is invalid. Resetting to 0, which will not filter based on max insert size.\n", config.maxIsize);
        config.maxIsize = 0;
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
