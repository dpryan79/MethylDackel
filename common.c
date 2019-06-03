#include "MethylDackel.h"
#include "version.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <limits.h>
#include <assert.h>
#include <pthread.h>

void parseBounds(char *s2, int *vals, int mult) {
    char *p, *s = strdup(s2), *end;
    int i, v;
    long tempV;

    p = strtok(s, ",");
    tempV = strtol(p, &end, 10);
    if((errno == ERANGE && (tempV == LONG_MAX || tempV == LONG_MIN)) || (errno != 0 && tempV == 0) || end == p) v = -1;
    else if(tempV > INT_MAX || tempV < LONG_MIN) v = -1;
    else v = tempV;

    if(v>=0) vals[4*mult] = v;
    else {
        fprintf(stderr, "Invalid bounds string, %s\n", s2);
        free(s);
        return;
    }
    for(i=1; i<4; i++) {
        p = strtok(NULL, ",");
        tempV = strtol(p, &end, 10);
        if((errno == ERANGE && (tempV == LONG_MAX || tempV == LONG_MIN)) || (errno != 0 && tempV == 0) || end == p) v = -1;
        else if(tempV > INT_MAX || tempV < LONG_MIN) v = -1;
        else v = tempV;

        if(v>=0) vals[4*mult+i] = v;
        else {
            fprintf(stderr, "Invalid bounds string, %s\n", s2);
            free(s);
            return;
        }
    }
    free(s);
}

void print_version() {
    printf("%s (using HTSlib version %s)\n", VERSION, hts_version());
}

inline int isCpG(char *seq, int pos, int seqlen) {
    if(*(seq+pos) == 'C' || *(seq+pos) == 'c') {
        if(pos+1 == seqlen) return 0;
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
    if(*(seq+pos) == 'C' || *(seq+pos) == 'c') {
        if(pos+2 >= seqlen) return 0;
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
    if(*(seq+pos) == 'C' || *(seq+pos) == 'c') return 1;
    else if(*(seq+pos) == 'G' || *(seq+pos) == 'g') return -1;
    return 0;
}

int getStrand(bam1_t *b) {
    char *XG = (char *) bam_aux_get(b, "XG");
    //Only bismark uses the XG tag like this. Some other aligners use it for other purposes...
    if(XG != NULL && *(XG+1) != 'C' && *(XG+1) != 'G') XG = NULL;
    if(XG == NULL) { //Can't handle non-directional libraries!
        if(b->core.flag & BAM_FPAIRED) {
            if((b->core.flag & 0x50) == 0x50) return 2; //Read1, reverse comp. == OB
            else if(b->core.flag & 0x40) return 1; //Read1, forward == OT
            else if((b->core.flag & 0x90) == 0x90) return 1; //Read2, reverse comp. == OT
            else if(b->core.flag & 0x80) return 2; //Read2, forward == OB
            return 0; //One of the above should be set!
        } else {
            if(b->core.flag & 0x10) return 2; //Reverse comp. == OB
            return 1; //OT
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

    if(strand==0) {
        fprintf(stderr, "Can't determine the strand of a read!\n");
        assert(strand != 0);
    }
    //Is the phred score even high enough?
    if(bam_get_qual(plp->b)[plp->qpos] < config->minPhred) return 0;

    if(base == 2 && (strand==1 || strand==3)) return 1; //C on an OT/CTOT alignment
    else if(base == 8 && (strand==1 || strand==3)) return -1; //T on an OT/CTOT alignment
    else if(base == 4 && (strand==2 || strand==4)) return 1; //G on an OB/CTOB alignment
    else if(base == 1 && (strand==2 || strand==4)) return -1; //A on an OB/CTOB alignment
    return 0;
}

//Convert bases outside of the bounds to N and their phred scores to 0
bam1_t *trimAlignment(bam1_t *b, int bounds[16]) {
    int strand = getStrand(b)-1;
    int i, lb, rb;
    uint8_t *qual = bam_get_qual(b);
    uint8_t *seq = bam_get_seq(b);

    if(b->core.flag & BAM_FREAD2) {
        lb = bounds[4*strand+2];
        rb = bounds[4*strand+3];
    } else {
        lb = bounds[4*strand];
        rb = bounds[4*strand+1];
    }

    lb = (lb<b->core.l_qseq) ? lb : b->core.l_qseq;

    //trim on the left
    if(lb) {
        for(i=0; i<lb; i++) {
            qual[i] = 0;
            if(i&1) seq[i>>1] |= 0xf;
            else seq[i>>1] |= 0xf0;
        }
    }

    //trim on the right
    if(rb) {
        for(i=rb; i<b->core.l_qseq; i++) {
            qual[i] = 0;
            if(i&1) seq[i>>1] |= 0xf;
            else seq[i>>1] |= 0xf0;
        }
    }

    return b;
}

bam1_t *trimAbsoluteAlignment(bam1_t *b, int bounds[16]) {
    int strand = getStrand(b)-1;
    int i, lb, rb;
    uint8_t *qual = bam_get_qual(b);
    uint8_t *seq = bam_get_seq(b);

    if(b->core.flag & BAM_FREAD2) {
        lb = bounds[4*strand+2];
        rb = bounds[4*strand+3];
    } else {
        lb = bounds[4*strand];
        rb = bounds[4*strand+1];
    }

    lb = (lb<b->core.l_qseq) ? lb : b->core.l_qseq;
    rb = (rb<b->core.l_qseq) ? rb : b->core.l_qseq;

    if(lb) {
        for(i=0; i<lb; i++) {
            qual[i] = 0;
            if(i&1) seq[i>>1] |= 0xf;
            else seq[i>>1] |= 0xf0;
        }
    }

    if(rb) {
        for(i=0; i<rb; i++) {
            qual[b->core.l_qseq - i] = 0;
            if((b->core.l_qseq - i)&1) seq[(b->core.l_qseq - i)>>1] |= 0xf;
            else seq[(b->core.l_qseq - i)>>1] |= 0xf0;
        }
    }

    return b;
}

char* getMappabilityValue(Config* config, char* chrom_n, uint32_t start, uint32_t end)
{
    fprintf(stderr, "started getMappabilityValue\n");
    uint32_t chrom = -1;
    for(int i = 0; i<config->BW_ptr->cl->nKeys; i++) //loop over chromosomes
    {
        //fprintf(stderr, "getting chrom id\n");
        if(!strcmp(config->BW_ptr->cl->chrom[i], chrom_n)) //found the chromosome
        {
            //fprintf(stderr, "found chrom id\n");
            chrom = i;
            break;
        }
    }
    fprintf(stderr, "out of chrom id loop\n");
    char* data = malloc((end-start)*sizeof(char)); //allocate array for data
    fprintf(stderr, "data array allocated\n");
    int index = start/8;
    int offset = start%8;
    fprintf(stderr, "calculated index and offset\n");
    for(int i = 0; i<end-start; i++)
    {
        //fprintf(stderr, "looping over data\n");
        char byte = config->bw_data[chrom][index];
        //fprintf(stderr, "got data byte\n");
        char mask = 1 >> offset;
        //fprintf(stderr, "created mask\n");
        char val = (byte & mask) << offset;
        //fprintf(stderr, "masked data byte and did offset\n");
        data[i] = val;
        //fprintf(stderr, "assigned to data arr\n");
        if(offset == 7)
        {
            //fprintf(stderr, "moving to next byte\n");
            index++;
            offset = 0;
        }
        else
        {
            //fprintf(stderr, "incrementing offset\n");
            offset++;
        }
        
    }
    fprintf(stderr, "done with getting data\n");
    return data;
}

char check_mappability(void *data, bam1_t *b) {
    //returns number of mappable reads in read pair (0-2)
    mplp_data *ldata = (mplp_data *) data;
    int read1_start;
    int read1_end;
    int read2_start;
    int read2_end;
    int i; //loop index
    int num_mappable_reads = 0;
    char num_mappable_bases = 0; //counter
    char *vals = NULL;
    if(b->core.flag & BAM_FREAD1 || (bam_is_rev(b) && b->core.flag & BAM_FREAD2)) //is this the left read?
    {
        read1_start = b->core.pos;
        read1_end = b->core.pos + b->core.l_qseq;
        read2_start = b->core.mpos;
        read2_end = b->core.mpos + b->core.l_qseq; //assuming both reads same length to avoid issues finding read2_end on right read (doing the same on the left read for consistency)
        
    }
    else //get pos for right read
    {
        read2_start = b->core.pos;
        read2_end = b->core.pos + b->core.l_qseq;
        read1_start = b->core.mpos;
        read1_end = b->core.mpos + b->core.l_qseq; //assuming both reads same length to avoid issues finding read2_end
    }
    if(ldata->config->BW_ptr == NULL) //invalid or missing bigWig
    {
        return -1; //this is checked in filter_func as well, so this statement should never run
    }
    //pthread_mutex_lock(&bwMutex); //locking to avoid threading issues on read
    vals = getMappabilityValue(ldata->config, ldata->hdr->target_name[b->core.tid], read1_start, read1_end+1);
    //bwGetValues(ldata->config->BW_ptr, ldata->hdr->target_name[b->core.tid], read1_start, read1_end+1, 1);
    //pthread_mutex_unlock(&bwMutex);
    
    for (i=0; i<=read1_end-read1_start; i++)
    {
        if(vals[i]) //considering NaN as 0 so as to not call reads mappable unless there is data saying they are
        {
            num_mappable_bases++;
        }
        if(num_mappable_bases >= ldata->config->minMappableBases)
        {
            num_mappable_reads++;
            break; //done with this read
        }
    }
    //free(vals->start);
    //free(vals->end);
    //free(vals->value);
    free(vals);
    //pthread_mutex_lock(&bwMutex); //locking to avoid threading issues on read
    vals = getMappabilityValue(ldata->config, ldata->hdr->target_name[b->core.tid], read2_start, read2_end+1);
    //bwGetValues(ldata->config->BW_ptr, ldata->hdr->target_name[b->core.tid], read2_start, read2_end+1, 1);
    //pthread_mutex_unlock(&bwMutex);
    
    num_mappable_bases = 0;
    for (i=0; i<=read2_end-read2_start; i++)
    {
        if(vals[i]) //considering NaN as 0
        {
            num_mappable_bases++;
        }
        if(num_mappable_bases >= ldata->config->minMappableBases)
        {
            num_mappable_reads++;
            break; //done with this read
        }
    }
    //free(vals->start);
    //free(vals->end);
    //free(vals->value);
    free(vals);
    return num_mappable_reads;
}

//This will need to be restructured to handle multiple input files
int filter_func(void *data, bam1_t *b) {
    int rv, NH, overlap;
    mplp_data *ldata = (mplp_data *) data;
    uint8_t *p;

    while(1) {
        rv = ldata->iter ? sam_itr_next(ldata->fp, ldata->iter, b) : sam_read1(ldata->fp, ldata->hdr, b);

        if(rv<0) return rv;
        if(b->core.tid == -1 || b->core.flag & BAM_FUNMAP) continue; //Unmapped
        if(b->core.qual < ldata->config->minMapq) continue; //-q
        if(b->core.flag & ldata->config->ignoreFlags) continue; //By default: secondary alignments, QC failed, PCR duplicates, and supplemental alignments
        if(ldata->config->requireFlags && (b->core.flag & ldata->config->requireFlags) != ldata->config->requireFlags) continue;
        if(!ldata->config->keepDupes && b->core.flag & BAM_FDUP) continue;
        p = bam_aux_get(b, "NH");
        if(p != NULL) {
            NH = bam_aux2i(p);
            if(NH>1) continue; //Ignore obvious multimappers
        }
        //check_mappability(ldata, b);
        //printf("%i", !ldata->config->BW_ptr);
        if(ldata->config->BW_ptr && check_mappability(ldata, b) == 0) continue; //Low mappability
        if(!ldata->config->keepSingleton && (b->core.flag & 0x9) == 0x9) continue; //Singleton
        if(!ldata->config->keepDiscordant && (b->core.flag & 0x3) == 0x1) continue; //Discordant
        if((b->core.flag & 0x9) == 0x1) b->core.flag |= 0x2; //Discordant pairs can cause double counts
        if(ldata->config->bed) { //Prefilter reads overlapping a BED file (N.B., strand independent).
            overlap = spanOverlapsBED(b->core.tid, b->core.pos, bam_endpos(b), ldata->config->bed, &(ldata->bedIdx));
            if(overlap == 0) continue;
            if(overlap < 0) {
                rv = -1;
                break;
            }
        }

        /***********************************************************************
        *
        * Deal with bounds inclusion (--OT, --OB, etc.)
        * If we don't do this now, then dealing with this after the overlap
        * detection will result in losing a lot of calls that we actually should
        * keep (i.e., if a call is in an overlapping region near the end of read
        * #1 and that region is excluded in that read then we lose the call).
        * The overlap detection will only decrement read #2's phred score by 20%
        * (instead of to 0) if there's a base mismatch and read #2 has the
        * higher phred score at that position.
        *
        ***********************************************************************/
        if(ldata->config->bounds) b = trimAlignment(b, ldata->config->bounds);
        if(ldata->config->absoluteBounds) b = trimAbsoluteAlignment(b, ldata->config->absoluteBounds);
        break;
    }
    return rv;
}

//Ensure that CpGs and CHGs are never split between threads. Move end positions to the right
void adjustBounds(Config *config, bam_hdr_t *hdr, faidx_t *fai, uint32_t *localTid, uint32_t *localPos, uint32_t *localEnd) {
    uint32_t start, end, tmp; //For faidx_fetch_seq, these are 0-based fully closed!!!
    int seqlen;
    char *seq;

    end = *localEnd + 1;
    if(*localEnd > 0) {
        start = *localEnd - 1;
    } else {
        start = 0;
    }
    seq = faidx_fetch_seq(fai, hdr->target_name[*localTid], start, end, &seqlen);
    if(seqlen > 1) {
        if(seqlen > 2 && (seq[0] & 0x5F) == 'C' && (seq[2] & 0x5F) == 'G') { //CHG
            *localEnd += 2;
        } else if((seq[1] & 0x5F) == 'G') { //Possible CpG or CHG
            *localEnd += 1;
        }
    }
    free(seq);

    //though unlikely to ever not be the case, ensure start < end;
    if(*localPos > *localEnd) {
        tmp = *localPos;
        *localPos = *localEnd;
        *localEnd = tmp;
    }
}
