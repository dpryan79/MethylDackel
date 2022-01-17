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
    if(pos >= seqlen) return 0;
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
    if(pos >= seqlen) return 0;
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
    if(pos >= seqlen) return 0;
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

//Calculate the query moves for a specified reference distance (from the start of a read) 
int rlen2qlen(int rlen, int n_cigar, uint32_t *CIGAR, bam1_t *b) {
    int readPos = 0;
    int mapPos = 0;
    int cigarPos = 0;
    int cigarOPType;

    /* remarks
    In certain cases, the 3' ends of one read may extend beyond the 5' end of the mate (ie cases where insert size is ambiguous)
    For the purpose of trimming:
        - when query rlen > read rlen, the full read qlen is returned
        - when rlen is negative, return 0
    */
    if(rlen > bam_cigar2rlen(n_cigar, CIGAR)) {
        readPos = bam_cigar2qlen(n_cigar, CIGAR);
        mapPos = rlen;
    } else if(rlen < 0) {
        mapPos = rlen;
        fprintf(stderr, "Warning: rlen is negative; 0 qlen returned. qname is [[ %s ]].\n", bam_get_qname(b));
    }

    while(mapPos < rlen && cigarPos < n_cigar) {
        cigarOPType = bam_cigar_type(CIGAR[cigarPos]);
        if(cigarOPType & 2) {
            if(cigarOPType & 1) { // consumes both query and ref; move both mapPos and readPos
                if((mapPos + bam_cigar_oplen(CIGAR[cigarPos])) >= rlen) {
                    readPos = readPos + rlen - mapPos; // move readPos to the location where mapPos is
                    break;
                } else {
                    mapPos += bam_cigar_oplen(CIGAR[cigarPos]);
                    readPos += bam_cigar_oplen(CIGAR[cigarPos++]);
                    continue;
                } 
            } else { // consumes only ref but not query
                mapPos += bam_cigar_oplen(CIGAR[cigarPos++]);
                if(mapPos >= rlen) {
                    break; // readPos does not have to be moved
                } else {
                    continue;
                }
            }
        } else if(cigarOPType & 1) { // consumes only query but not ref
            readPos += bam_cigar_oplen(CIGAR[cigarPos++]);
            continue;
        } else { // consumes neither query or ref
            cigarPos++;
            continue;
        }
    }
    return(readPos);
}


//Calculate the least required query moves for a specified reference distance (from the start of a read) 
int qlen2rlen(int qlen, int n_cigar, uint32_t *CIGAR, bam1_t *b) {
    int readPos = 0;
    int mapPos = 0;
    int cigarPos = 0;
    int cigarOPType;

    /* remarks
    In the case where the desired trimming is greater than the read, return the full read rlen
    qlen should not be negative, but in that case, 0 is returned with a message to stderr
    */
    if(qlen > bam_cigar2qlen(n_cigar, CIGAR)) {
        mapPos = bam_cigar2rlen(n_cigar, CIGAR);
        readPos = qlen;
    } else if(qlen < 0) {
        readPos = qlen;
        fprintf(stderr, "Warning: qlen is negative; 0 rlen returned. qname is [[ %s ]].\n", bam_get_qname(b));
    }

    while(readPos < qlen && cigarPos < n_cigar) {
        cigarOPType = bam_cigar_type(CIGAR[cigarPos]);
        if(cigarOPType & 2) {
            if(cigarOPType & 1) { // consumes both query and ref; move both readPos and mapPos
                if((readPos + bam_cigar_oplen(CIGAR[cigarPos])) >= qlen) {
                    mapPos = mapPos + qlen - readPos; // move mapPos to the location where readPos is
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
            if(readPos >= qlen) {
                break; // mapPos does not have to be moved
            } else {
                continue;
            }
        } else { // consumes neither query or ref
            cigarPos++;
            continue;
        }
    }
    return(mapPos);
}

/* pseudocode
aim: to find out whether the other end needs trimming

approach:
calculate selfEndPos / mateEndPos;
do the two reads overlap?
    y:  calculate the mapPos of trimmed mate;
        does the trimming cross into current read?
            y:  trim for readLen required to reach mapPos
            n:  no need trim
    n:  is trimLen > mate qlen?
            y:  tricky case 1, need to assume that the observed isize is the true fragment length, ie need to assume there is no indel in the uncovered area;
                calculate the hypothesized trim mapPos, and trim accordingly
            n:  no need trim

the reason to translate trim to mapPos first, then calculate qlen from startpos is that:
    mapPos is the only anchor between the paired reads (unless they overlap);
    and CIGAR string always go from left to right, and therefore must start from startpos

tricky case 2:
what's the interpretation when read is soft trimmed (ie bam_cigar_op(cigar) == BAM_CSOFT_CLIP)?
    custom barcode? adapter?
    at this point, I'll stop worrying about them and just use the read-5' ends as the gold standard
*/


bam1_t *trimFragmentEnds(bam1_t *b, int fivePrime, int threePrime) {
    int i, lb = 0, rb = 0, selfTrim, mateTrim;
    int mateTrim_mapPos;
    uint8_t *qual = bam_get_qual(b);
    uint8_t *seq = bam_get_seq(b);

    // parsing mate cigar string
    char *m_cstring = bam_aux2Z(bam_aux_get(b, "MC"));
    char *end;
    uint32_t *CIGAR = bam_get_cigar(b);
    uint32_t *m_CIGAR = NULL;
    size_t m = 0;
    int m_n_cigar;
    m_n_cigar = sam_parse_cigar(m_cstring, &end, &m_CIGAR, &m);

    uint32_t mateEndPos = b->core.mpos + bam_cigar2rlen(m_n_cigar, m_CIGAR);
    uint32_t selfEndPos = bam_endpos(b);
    int m_qlen = bam_cigar2qlen(m_n_cigar, m_CIGAR);

    // set the params
    if(b->core.flag & BAM_FPROPER_PAIR) {
        selfTrim = (b->core.flag & BAM_FREAD1) ? fivePrime : threePrime;
        mateTrim = (b->core.flag & BAM_FREAD1) ? threePrime : fivePrime;
        
        if(b->core.flag & BAM_FMREVERSE) { // OT; selfTrim at left bound, mateTrim at right bound
            lb = selfTrim;
            if(selfEndPos > b->core.mpos) { // ie reads overlap
                if(m_qlen >= mateTrim) {
                    mateTrim_mapPos = b->core.mpos + qlen2rlen(m_qlen - mateTrim, m_n_cigar, m_CIGAR, b);
                    rb = (selfEndPos > mateTrim_mapPos) ? b->core.l_qseq - rlen2qlen(mateTrim_mapPos - b->core.pos, b->core.n_cigar, CIGAR, b) : 0;
                } else rb = b->core.l_qseq - rlen2qlen(b->core.mpos - b->core.pos, b->core.n_cigar, CIGAR, b) + (mateTrim - m_qlen);
            } else { // ie not overlap
                if(mateTrim > m_qlen) { // ie tricky case 1
                    mateTrim_mapPos = b->core.mpos - (mateTrim - m_qlen);
                    rb = (selfEndPos > mateTrim_mapPos) ? b->core.l_qseq - rlen2qlen(mateTrim_mapPos - b->core.pos, b->core.n_cigar, CIGAR, b) : 0;
                } else rb = 0;
            }
        } else { // OB; selfTrim at right bound, mateTrim at left bound
            rb = selfTrim;
            if(mateEndPos > b->core.pos) { // ie reads overlap
                if(m_qlen >= mateTrim) {
                    mateTrim_mapPos = b->core.mpos + qlen2rlen(mateTrim, m_n_cigar, m_CIGAR, b);
                    lb = (mateTrim_mapPos > selfEndPos) ? rlen2qlen(mateTrim_mapPos - b->core.pos, b->core.n_cigar, CIGAR, b) : 0;
                } else lb = rlen2qlen(mateEndPos - b->core.pos, b->core.n_cigar, CIGAR, b) + (mateTrim - m_qlen);
            } else { // ie not overlap
                if(mateTrim > m_qlen) { // ie tricky case 1
                    mateTrim_mapPos = mateEndPos + (mateTrim - m_qlen);
                    lb = (mateTrim_mapPos > selfEndPos) ? rlen2qlen(mateTrim_mapPos - b->core.pos, b->core.n_cigar, CIGAR, b) : 0;
                } else lb = 0;
            }
        }
    } else {
        fprintf(stderr, "Read with qname [[ %s ]] not properly paired\n", bam_get_qname(b));
    }

    lb = (lb < b->core.l_qseq) ? lb : b->core.l_qseq;
    rb = (rb < b->core.l_qseq) ? rb : b->core.l_qseq;

    if(lb) {
        for(i=0; i<lb; i++) {
            qual[i] = 0;
            if(i&1) seq[i>>1] |= 0xf; // 0xf takes 4 bit, corresponds to N
            else seq[i>>1] |= 0xf0; // 1 byte of seq encode 2 base, and 1 byte of qual encode 1 base qual, thus the shift is needed
        }
    }

    if(rb) {
        for(i=b->core.l_qseq - rb; i<b->core.l_qseq; i++) {
            qual[i] = 0;
            if(i&1) seq[i>>1] |= 0xf;
            else seq[i>>1] |= 0xf0;
        }
    }
    free(m_CIGAR);
    return b;
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

unsigned char* getMappabilityValue(Config* config, char* chrom_n, uint32_t start, uint32_t end)
{
    char chromFound = 0;
    uint32_t chrom = -1;
    int i;
    for(i = 0; i<config->chromCount; i++) //loop over chromosomes
    {
        if(!strcmp(config->chromNames[i], chrom_n)) //found the chromosome
        {
            chrom = i;
            chromFound = 1;
            break;
        }
    }
    unsigned char* data = malloc((end-start)*sizeof(unsigned char)); //allocate array for data
    int index = start/8;
    int offset = start%8;
    //int startindex = start/8; //debug
    //int startoffset = start%8; //debug

    //debug
    int arrlen = 0; //variable to store the length of the array used for the data (this is not the same as the chromosome length, as each value is one bit and this is an array of characters, i.e. bytes)
    if(chromFound)
    {
        arrlen = config->chromLengths[chrom]/8; //array length is chromosome length over 8 (number of bits to number of bytes)
        if(config->chromLengths[chrom]%8 > 0) //if there is a remainder that didn't divide evenly
        {
            arrlen++; //add an extra byte to store it
        }
    }

    for(i = 0; i<end-start; i++)
    {
        unsigned char byte;
        if(chromFound) //was a chrom ID found for chrom_n, or is chrom still -1 (or here 4,294,967,295) i.e. chrom not found
        {
            if(index>=arrlen)
            {
                byte = 0;
                //printf("warning: accessing at index %d offset %d, arrlen is %d!\ni is %d, end-start is %d, end pos should be %d!\nstartindex is %d, startoffset is %d, start is %d, end is %d\nchrom len is %d\n", index, offset, arrlen, i, end-start, (start/8)+((end-start)/8), startindex, startoffset, start, end, config->chromLengths[chrom]);
            }
            else
            {
                byte = config->bw_data[chrom][index];
            }
        }
        else
        {
            byte = 0; //if not a valid chrom, mappability is N/A i.e. 0
        }
        unsigned char mask = 1 << offset;
        unsigned char val = (byte & mask) >> offset;
        data[i] = val;
        if(offset == 7) //at end of byte
        {
            index++;
            offset = 0;
        }
        else
        {
            offset++;
        }
        
    }
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
    unsigned char *vals = NULL;
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
    vals = getMappabilityValue(ldata->config, ldata->hdr->target_name[b->core.tid], read1_start, read1_end);
    
    for (i=0; i<read1_end-read1_start; i++)
    {
        if(vals[i] > 0) //is base above threshold?
        {
            num_mappable_bases++;
        }
        if(num_mappable_bases >= ldata->config->minMappableBases)
        {
            num_mappable_reads++;
            break; //done with this read
        }
    }
    free(vals);
    vals = getMappabilityValue(ldata->config, ldata->hdr->target_name[b->core.tid], read2_start, read2_end);
    
    num_mappable_bases = 0;
    for (i=0; i<read2_end-read2_start; i++)
    {
        if(vals[i] > 0) //is base above threshold?
        {
            num_mappable_bases++;
        }
        if(num_mappable_bases >= ldata->config->minMappableBases)
        {
            num_mappable_reads++;
            break; //done with this read
        }
    }
    free(vals);
    return num_mappable_reads;
}

// This is the same as updateMetrics, 1 on methylation, -1 on unmethylation
int getMethylState(bam1_t *b, int seqPos, Config *config) {
    uint8_t base = bam_seqi(bam_get_seq(b), seqPos);
    int strand = getStrand(b); //1=OT, 2=OB, 3=CTOT, 4=CTOB

    if(strand==0) {
        fprintf(stderr, "Can't determine the strand of a read!\n");
        assert(strand != 0);
    }
    //Is the phred score even high enough?
    if(bam_get_qual(b)[seqPos] < config->minPhred) return 0;

    if(base == 2 && (strand==1 || strand==3)) return 1; //C on an OT/CTOT alignment
    else if(base == 8 && (strand==1 || strand==3)) return -1; //T on an OT/CTOT alignment
    else if(base == 4 && (strand==2 || strand==4)) return 1; //G on an OB/CTOB alignment
    else if(base == 1 && (strand==2 || strand==4)) return -1; //A on an OB/CTOB alignment
    return 0;
}

float computeEfficiency(unsigned int nMethyl, unsigned int nUMethyl) {
    if(nMethyl + nUMethyl == 0) return 1.0;
    return nUMethyl / ((float)(nMethyl + nUMethyl));
}

float computeConversionEfficiency(bam1_t *b, mplp_data *ldata) {
    unsigned int nMethyl = 0, nUMethyl = 0;
    uint32_t i, j, seqEnd = ldata->offset + ldata->lseq;  // 1-base after the end of the sequence
    uint32_t *cigar = bam_get_cigar(b), op, opLen;
    int direction, state;
    int pos = b->core.pos, seqPos = 0;  //position in the genome and position in the read

    for(i=0; i<b->core.n_cigar; i++) {
        op = bam_cigar_op(cigar[i]);
        opLen = bam_cigar_oplen(cigar[i]);

        switch(op) {
        case 0:
        case 7:
        case 8:
            // do something
            for(j=0; j<opLen; j++, seqPos++) {
                if(pos+j >= seqEnd) return computeEfficiency(nMethyl, nUMethyl);
                if(isCpG(ldata->seq, pos+j-ldata->offset, ldata->lseq)) {
                    continue;
                } else if((direction = isCHG(ldata->seq, pos+j-ldata->offset, ldata->lseq))) {
                    state = getMethylState(b, seqPos, ldata->config);
                    if(state > 0) nMethyl++;
                    else if(state < 0) nUMethyl++;
                } else if((direction = isCHH(ldata->seq, pos+j-ldata->offset, ldata->lseq))) {
                    state = getMethylState(b, seqPos, ldata->config);
                    if(state > 0) nMethyl++;
                    else if(state < 0) nUMethyl++;
                }
            }
            break;
        case 1:  // I, consume read
        case 4:  // S, consume read
            seqPos += opLen;
            break;
        case 2:  // D, consume seq
        case 3:  // N, consume seq
            pos += opLen;
            break;
        }
    }

    return computeEfficiency(nMethyl, nUMethyl);
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
        if(ldata->config->minIsize && abs(b->core.isize) < ldata->config->minIsize) continue; //Minimum insert size
        if(ldata->config->maxIsize && abs(b->core.isize) > ldata->config->maxIsize) continue; //Maximum insert size
        if(b->core.flag & ldata->config->ignoreFlags) continue; //By default: secondary alignments, QC failed, PCR duplicates, and supplemental alignments
        if(ldata->config->requireFlags && (b->core.flag & ldata->config->requireFlags) != ldata->config->requireFlags) continue;
        if(!ldata->config->keepDupes && b->core.flag & BAM_FDUP) continue;
        if(!ldata->config->ignoreNH) {
            p = bam_aux_get(b, "NH");
            if(p != NULL) {
                NH = bam_aux2i(p);
                if(NH>1) continue; //Ignore obvious multimappers
            }
        }
        if((ldata->config->filterMappability) && check_mappability(ldata, b) == 0) continue; //Low mappability
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

        // This is (A) moderately expensive to compute and (B) not completely correct at chunk boundaries.
        if(ldata->config->minConversionEfficiency > 0.0) {
            if(computeConversionEfficiency(b, ldata) < ldata->config->minConversionEfficiency) continue;
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
        if(ldata->config->fivePrime || ldata->config->threePrime) b = trimFragmentEnds(b, ldata->config->fivePrime, ldata->config->threePrime);
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
