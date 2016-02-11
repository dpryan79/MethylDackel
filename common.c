#include "PileOMeth.h"
#include "version.h"
#include <assert.h>

void print_version() {
    printf("%s (using HTSlib version %s)\n", VERSION, hts_version());
}

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

        /***********************************************************************
        *
        * Deal with bounds inclusion (--OT, --OB, etc.)
        * If we don't do this now, then the dealing with this after the overlap
        * detection will result in losing a lot of calls that we actually should
        * keep (i.e., if a call is in an overlapping region near an end of read
        * #1 and that region is excluded in that read then we lose the call).
        * The overlap detection will only decrement read #2's phred score by 20%
        * (instead of to 0) if there's a base mismatch and read #2 has the
        * higher phred score at that position.
        *
        ***********************************************************************/
        if(ldata->config->bounds) b = trimAlignment(b, ldata->config->bounds);
        break;
    }
    return rv;
}
