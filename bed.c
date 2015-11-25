#include "PileOMeth.h"
#include <ctype.h>
#include <limits.h>
#include "htslib/kseq.h"
#include <assert.h>
KSTREAM_INIT(gzFile, gzread, 8192)

//returns <0 if region0 comes before region1
//returns >0 if region0 comes after region1
//returns 0 on overlap
inline int64_t compareRegions(int32_t tid0, int32_t start0, int32_t end0, int32_t tid1, int32_t start1, int32_t end1) {
    if(tid0 != tid1) return ((int64_t) tid0)-((int64_t) tid1);
    if(start0 < start1 && end0 >= start1) return 0;
    if(start0 >= start1 && start0 < end1) return 0;
    return ((int64_t) start0)-((int64_t) start1);
}
    
//1: Overlaps
//0: No Overlap
//-1: End of BED file
//This should eventually be modified to use a hash
int spanOverlapsBED(int32_t tid, int32_t start, int32_t end, bedRegions *regs, int *idx) {
    bedRegion *reg = regs->region;
    int64_t rv = 2;
    int i;

    //First, test the last position we probed.
    if(compareRegions(reg[*idx].tid, reg[*idx].start, reg[*idx].end-1, tid, start, end) == 0) return 1;
    else {
        for(i=*idx; i<regs->n; i++) {
            rv = compareRegions(reg[i].tid, reg[i].start, reg[i].end-1, tid, start, end);
            if(rv >= 0) {
                *idx = i;
                rv = (rv >= 1) ? 0 : 1;
                break;
            }
        }
        if(rv < 0) rv = -1;
    }
    assert(rv!=2);
    return rv;
}

//-1 if region is before pos (need to go to the next region)
// 0 if pos not in region
//>0 if pos in region
int posOverlapsBED(int32_t tid, int32_t pos, bedRegions *regions, int idx) {
    if(idx >= regions->n) return 0; //Past the end of the regions

    if(tid != regions->region[idx].tid) return (regions->region[idx].tid < tid) ? -1 : 0;
    if(pos >= regions->region[idx].end) return -1;
    if(pos < regions->region[idx].start) return 0;
    return 1;
}

//Return 1 on overlap, otherwise 0
int readStrandOverlapsBED(bam1_t *b, bedRegion region) {
    int s = getStrand(b);
    if(region.strand) {
        if(region.strand == 1 && (s==1 || s==3)) return 1;
        if(region.strand == 2 && (s==2 || s==4)) return 1;
        return 0;
    }
    return 1;
}

int sortBED_func(const void *a, const void *b) {
    bedRegion *pa = (bedRegion *) a;
    bedRegion *pb = (bedRegion *) b;

    if(pa->tid < pb->tid) return -1;
    if(pa->tid > pb->tid) return 1;
    if(pa->start < pb->start) return -1;
    if(pa->start > pb->start) return 1;
    if(pa->end < pb->end) return -1;
    if(pa->end > pb->end) return 1;
    //This is arbitrary
    if(pa->strand < pb->strand) return -1;
    if(pa->strand > pb->strand) return 1;
    return 0; //This shouldn't actually happen
}

//For convenience
void sortBED(bedRegions *regions) {
    qsort((void *) regions->region, regions->n, sizeof(bedRegion), sortBED_func);
}

//Returns a NULL pointer on error, which should cause the program to exit
//N.B., If column 4 exists it can't contain a space (quotes won't help)
//It might make sense to just use a regex
bedRegions *parseBED(char *fn, bam_hdr_t *hdr) {
    gzFile fp = NULL;
    int i, dret = 0; //dret isn't actually used
    char *p1, *p2;
    int32_t lnum = 0;
    kstream_t *ks = NULL;
    kstring_t str = {0, 0, NULL};
    bedRegions *regions = NULL;

    if((fp = gzopen(fn, "r")) == NULL) {
        fprintf(stderr, "Couldn't open %s for reading.\n", fn);
        goto err;
    }
    if((regions = malloc(sizeof(bedRegions))) == NULL) {
        fprintf(stderr, "Couldn't allocate sufficient space to hold the BED file.\n");
        goto err;
    }
    ks = ks_init(fp);

    //Initialize the regions struct
    regions->n = 0;
    regions->m = 1000;
    if((regions->region = malloc(regions->m * sizeof(bedRegion))) == NULL) {
        fprintf(stderr, "Couldn't allocate enough space to store regions in %s.\n Too little memory?", fn);
        goto err;
    }

    //Iterate over the lines, skipping apparent comments and header lines
    while(ks_getuntil(ks, KS_SEP_LINE, &str, &dret) > 0) {
        lnum++;
        p1 = str.s;
        p2 = p1;

        //Skip blank lines, comments and any line beginning with "track" or "browser"
        if(*p1 == '\0') continue;
        if(*p1 == '#') continue;

        //Do we need to grow regions?
        if(regions->m - regions->n < 100) {
            regions->m += 1000;
            if((regions->region = realloc(regions->region, regions->m * sizeof(bedRegion))) == NULL) {
                fprintf(stderr, "We couldn't increase the space occupied by the regions structure to hold %" PRId32 " elements!\n", regions->m);
                fprintf(stderr, "You might need more memory\n");
                goto err;
            }
        }

        //Initialize the region
        regions->region[regions->n].tid = -1;
        regions->region[regions->n].start = -1;
        regions->region[regions->n].end = -1;
        regions->region[regions->n].strand = 0;

        //Parse the chromosome/contig name
        while(*p2 && !isspace(*p2)) p2++;
        if(*p2 != '\0') *p2 = '\0';
        for(i=0; i<hdr->n_targets; i++) {
            if(strcmp(p1, hdr->target_name[i]) == 0) {
                regions->region[regions->n].tid = i;
                break;
            }
        }
        //Just in case someone thought calling a contig "track" or "browser" was a good idea...
        if(regions->region[regions->n].tid == -1) {
            if(strcmp(p1, "track") == 0) continue;
            if(strcmp(p1, "browser") == 0) continue;
        }
        /***********************************************************************
        *
        * If we couldn't parse the region then we'll count this as an error.
        * It's possible that there's something we overlooked, since the BED
        * format is poorly specified. Consequently, future versions might simply
        * skip such lines.
        *
        ***********************************************************************/
        if(regions->region[regions->n].tid == -1) {
            fprintf(stderr, "Couldn't properly parse line number %i in %s.\n", lnum, fn);
            goto err;
        }

        //Parse the start and end positions
        p1 = p2+1;
        if(sscanf(p1, "%"PRId32, &regions->region[regions->n].start) != 1 || regions->region[regions->n].start == -1) {
            fprintf(stderr, "Line %" PRId32 " of %s is malformed.\n", lnum, fn);
            goto err;
        }
        p2++;

        while(*p2 && !isspace(*p2)) p2++;
        if(*p2 != '\0') *p2 = '\0';
        p1 = p2+1;
        if(sscanf(p1, "%"PRId32, &regions->region[regions->n].end) != 1 || regions->region[regions->n].end == -1) {
            fprintf(stderr, "Line %" PRId32 " of %s is malformed.\n", lnum, fn);
            goto err;
        }
        if(regions->region[regions->n].start >= regions->region[regions->n].end) {
            fprintf(stderr, "The position on line %" PRId32 " of %s is incorrect (%"PRId32" >= %"PRId32".\n", lnum, fn, regions->region[regions->n].start, regions->region[regions->n].end);
            goto err;
        }
        if(regions->region[regions->n].start < 0) regions->region[regions->n].start = 0;
        if(regions->region[regions->n].end > hdr->target_len[regions->region[regions->n].tid]+1) \
            regions->region[regions->n].end = hdr->target_len[regions->region[regions->n].tid]+1;

        //It's possible to have an additional 3 (or more) fields, in which case we need to skip 2 and parse the third
        regions->n++;
        if(!isblank(*p2)) {
            continue;
        }

        //4th column
        while(*p2 && isspace(*p2)) p2++; //skip white-space
        if(*p2 == '\0') continue; //Just some excess white-space
        while(*p2 && !isspace(*p2)) p2++;
        if(*p2 == '\0') continue;

        //5th column
        while(*p2 && isspace(*p2)) p2++;
        if(*p2 == '\0') continue;
        while(*p2 && !isspace(*p2)) p2++;
        if(*p2 == '\0') continue;

        //Strand
        while(*p2 && isspace(*p2)) p2++;
        if(*p2 == '\0') continue;
        if(*p2 == '+') regions->region[regions->n-1].strand = 1;
        else if(*p2 == '-') regions->region[regions->n-1].strand = 2;
        continue;
    }

    ks_destroy(ks);
    free(str.s);
    gzclose(fp);

    //In case the order isn't as expected...
    sortBED(regions);

    fprintf(stderr, "Parsed %" PRId32 " regions in %s\n", regions->n, fn);
    return regions;

err :
    if(regions) destroyBED(regions);
    if(ks) ks_destroy(ks);
    if(str.s) free(str.s);
    if(fp) gzclose(fp);
    return NULL;
}

void destroyBED(bedRegions *regions) {
    free(regions->region);
    free(regions);
}
