#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>
#include <inttypes.h>
#include "htslib/faidx.h"
#include "htslib/kseq.h"

void print_version(void);

int fgets2(FILE*f, unsigned char *s, int l) {
    return (fgets((char*) s, l, f))?strlen((char*)s):0;
}

KSTREAM_INIT(FILE*, fgets2, 65536)

struct lastCall{
    char *chrom;
    int32_t start, end;
    uint32_t nmethyl, nunmethyl;
};

inline void printRecord(FILE *of, char *chr, int32_t start, int32_t end, uint32_t nmethyl, uint32_t nunmethyl) {
    fprintf(of, "%s\t%"PRId32"\t%"PRId32"\t%i\t%"PRIu32"\t%"PRIu32"\n", chr, \
        start, end, (int) (100.0 * ((double) nmethyl)/(nmethyl+nunmethyl)), \
        nmethyl, nunmethyl);
}

void MergeOrPrint(FILE *of, struct lastCall *last, char *chr, int32_t start, int32_t width, uint32_t nmethyl, uint32_t nunmethyl) {
    int32_t end;

    if(width>0) {
        end = start+width;
    } else {
        end = start+1;
        start = end+width;
    }

    if(last->chrom && strcmp(last->chrom, chr)==0 && last->start==start && last->end==end) {
        printRecord(of, chr, start, end, nmethyl+last->nmethyl, nunmethyl+last->nunmethyl);
        free(last->chrom);
        free(chr);
        last->chrom=NULL;
    } else {
        if(last->chrom) {
            printRecord(of, last->chrom, last->start, last->end, last->nmethyl, last->nunmethyl);
            free(last->chrom);
        }
        last->chrom = chr;
        last->start = start;
        last->end = end;
        last->nmethyl = nmethyl;
        last->nunmethyl = nunmethyl;
    }
}

int getContext(faidx_t *fai, char *chr, int32_t pos, int *width) {
    int len = faidx_seq_len(fai, chr);
    int rv = 2;
    int32_t start = (pos>2) ? pos-2 : 0;
    int32_t end = (pos+2<len) ? pos+2 : len-1;
    char *seq = faidx_fetch_seq(fai, chr, start, end, &len);
    if(!seq) return 3;
    char *base = seq+(pos-start);

    if(toupper(*base) == 'C') {
        if(end-pos) {
            if(toupper(*(base+1)) == 'G') {
                *width = 2;
                rv = 0;
            } else if(end-pos==2) {
                if(toupper(*(base+2)) == 'G') {
                    *width = 3;
                    rv = 1;
                }
            }
        }
    } else {
        assert(toupper(*base) == 'G');
        if(pos-start) {
            if(toupper(*(base-1)) == 'C') {
                *width = -2;
                rv = 0;
            } else if(pos-start==2) {
                if(toupper(*(base-2)) == 'C') {
                    *width = -3;
                    rv = 1;
                }
            }
        }
    }

    free(seq);
    return rv;
}

void mergeContext(FILE *ifile, faidx_t *fai, FILE *ofile) {
    struct lastCall *lastCpG = calloc(1, sizeof(struct lastCall));
    struct lastCall *lastCHG = calloc(1, sizeof(struct lastCall));
    assert(lastCpG);
    assert(lastCHG);
    char *p, *p2;
    int ret, type, width;
    char *chr;
    int32_t start, end;
    uint32_t nmethyl, nunmethyl;
    kstring_t str;
    str.s = 0; str.l = str.m = 0;
    kstream_t *ks = ks_init(ifile);
    assert(ks);

    while(ks_getuntil(ks, KS_SEP_LINE, &str, &ret) >= 0) {
        assert(str.l>0);
        if(strncmp(str.s, "track", 5) == 0) continue;
        p = strtok(str.s, "\t");
        chr = strdup(p);
        p = strtok(NULL, "\t");
        start = strtoll(p, &p2, 10);
        assert(p2!=p);
        p = strtok(NULL, "\t");
        end = strtoll(p, &p2, 10);
        assert(p2!=p);
        p = strtok(NULL, "\t");
        p = strtok(NULL, "\t");
        nmethyl = strtoul(p, &p2, 10);
        assert(p2!=p);
        p = strtok(NULL, "\n");
        nunmethyl = strtoul(p, &p2, 10);
        assert(p2!=p);

        type = getContext(fai, chr, start, &width);
        if(type==0) {
            MergeOrPrint(ofile, lastCpG, chr, start, width, nmethyl, nunmethyl);
        } else if(type==1) {
            MergeOrPrint(ofile, lastCHG, chr, start, width, nmethyl, nunmethyl);
        } else if(type==2) {
            printRecord(ofile, chr, start, end, nmethyl, nunmethyl);
            free(chr);
        } else {
            fprintf(stderr, "[mergeContext] Error, %s is an unknown chromosome name!\n", chr);
            free(chr);
            break;
        }
    }

    if(lastCpG->chrom) {
        printRecord(ofile, lastCpG->chrom, lastCpG->start, lastCpG->end, lastCpG->nmethyl, lastCpG->nunmethyl);
        free(lastCpG->chrom);
    }
    if(lastCHG->chrom) {
        printRecord(ofile, lastCHG->chrom, lastCHG->start, lastCHG->end, lastCHG->nmethyl, lastCHG->nunmethyl);
        free(lastCHG->chrom);
    }
    if(str.m) free(str.s);
    ks_destroy(ks);
    free(lastCpG);
    free(lastCHG);
}
void mergeContext_usage() {
    fprintf(stderr, "\nUsage: PileOMeth mergeContext [OPTIONS] <ref.fa> <input>\n");
    fprintf(stderr,
"\n"
"This program will merge single Cytosine methylation metrics into per-CpG/CHG\n"
"metrics. The input file must be coordinate sorted but is not required to contain\n"
"only a single sequence context, though not doing so may result in an unsorted\n"
"result.\n"
"\n"
"  ref.fa    Reference genome in fasta format. This must be indexed with\n"
"            samtools faidx\n"
"  input     An input file such as that produced by PileOMeth extract. Specifying\n"
"            - allows reading from a pipe.\n"
"\nOptions:\n"
"  -o STR    Output file name [stdout]\n"
"  --version Printer version and quit\n"
);
}

int mergeContext_main(int argc, char *argv[]) {
    faidx_t *fai;
    FILE *ifile, *ofile = stdout;
    char c;

    static struct option lopts[] = {
        {"help",    0, NULL, 'h'},
        {"version", 0, NULL, 'v'},
        {0,         0, NULL,   0}
    };
    while((c = getopt_long(argc, argv, "hvo:", lopts, NULL)) >= 0) {
        switch(c) {
        case 'h' :
            mergeContext_usage();
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
        default :
            fprintf(stderr, "Invalid option '%c'\n", c);
            mergeContext_usage();
            return 1;
        }
    }

    if(argc == 1) {
        mergeContext_usage();
        return 0;
    }
    if(argc - optind != 2) {
        fprintf(stderr, "You must supply a reference genome in fasta format and an input bedGraph files\n");
        mergeContext_usage();
        return -1;
    }
    if((fai = fai_load(argv[optind])) == NULL) {
        fprintf(stderr, "Couldn't open the index for %s!\n", argv[optind]);
        mergeContext_usage();
        return -2;
    }
    if((ifile = fopen(argv[optind+1], "r")) == NULL) {
        fprintf(stderr, "Couldn't open %s for reading!\n", argv[optind+1]);
        return -3;
    }

    fprintf(ofile, "track type=\"bedGraph\" description=\"merged Methylation metrics\"\n");
    mergeContext(ifile, fai, ofile);

    if(ofile != stdout) fclose(ofile);
    fclose(ifile);
    fai_destroy(fai);

    return 0;
}
