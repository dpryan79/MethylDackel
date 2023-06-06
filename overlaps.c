#include <inttypes.h>
#include <ctype.h>
#include <stdio.h>
#include <assert.h>
#include "htslib/khash.h"
#include "htslib/sam.h"
#include "MethylDackel.h"

int getStrand(bam1_t *b);

// Dictionary of overlapping reads
KHASH_MAP_INIT_STR(olap_hash, bam1_t *)
typedef khash_t(olap_hash) olap_hash_t;
khash_t(olap_hash) *ohash = NULL;


void * initOlapHash() {
    return (void*) kh_init(olap_hash);
}

void destroyOlapHash(void *ohash) {
    khash_t(olap_hash) * oh = ( khash_t(olap_hash) *) ohash;
    kh_destroy(olap_hash, oh);
}

//This is from bison
int32_t *calculate_positions(bam1_t *read) {
    int32_t *positions = malloc(sizeof(int32_t) * (size_t)read->core.l_qseq);
    int i, j, offset = 0, op, op_len;
    uint32_t *CIGAR = bam_get_cigar(read);
    int32_t previous_position = read->core.pos;

    for(i=0; i<read->core.n_cigar; i++) {
        op = bam_cigar_op(*(CIGAR+i));
        op_len = bam_cigar_oplen(*(CIGAR+i));
        for(j=0; j<op_len; j++) {
            if(op == 0 || op == 7 || op == 8) { //M, =, X
                *(positions+offset) = previous_position++;
                offset++;
            } else if(op == 1 || op == 4) { //I, S, H
                *(positions+offset) = -1;
                offset++;
            } else if(op == 2 || op == 3) { //D, N
                previous_position++;
            } else if(op == 5) { //H, which isn't in the sequence
            } else { //P
                fprintf(stderr, "[calculate_positions] We encountered a CIGAR operation that we're not ready to deal with in %s\n", bam_get_qname(read));
            }
        }
    }
    return positions;
}

static void cust_tweak_overlap_quality(bam1_t *a, bam1_t *b) {
    int ia = 0, ib = 0;
    int32_t na = a->core.l_qseq, nb = b->core.l_qseq;
    int32_t *posa = calculate_positions(a);
    int32_t *posb = calculate_positions(b);
    uint8_t *a_qual = bam_get_qual(a), *b_qual = bam_get_qual(b);
    uint8_t *a_seq  = bam_get_seq(a), *b_seq = bam_get_seq(b);

    //If alignments are on opposite strands then exit
    int sa = getStrand(a);
    int sb = getStrand(b);
    if(((sa-sb)&1) == 1) goto quit;

    //Go to the first mapped position
    while(ia<na && posa[ia]<0) ia++;
    while(ib<nb && posb[ib]<0) ib++;
    if(ia==na || ib==nb) goto quit;

    //Go to the first overlapping position
    if(posa[ia]<posb[ib]) {
        while(ia<na && posa[ia]<posb[ib]) ia++;
    } else {
        while(ib<nb && posb[ib]<posa[ia]) ib++;
    }
    if(ia==na || ib==nb) goto quit;

    //Take care of the overlap
    while(ia<na && ib<nb) {
        if(posa[ia] < posb[ib] || posa[ia] < 0) {
            ia++;
            continue;
        }
        if(posb[ib] < posa[ia] || posb[ib] < 0) {
            ib++;
            continue;
        }
        if(bam_seqi(a_seq, ia) != bam_seqi(b_seq, ib)) {
            if(a_qual[ia]>b_qual[ib] && bam_seqi(a_seq,ia) != 15) {
                a_qual[ia] -= b_qual[ib];
                b_qual[ib] = 0;
            } else if(b_qual[ib]>a_qual[ia] && bam_seqi(b_seq,ib) != 15) {
                b_qual[ib] -= a_qual[ia];
                a_qual[ia] = 0;
            } else {
                a_qual[ia] = 0;
                b_qual[ib] = 0;
            }
        } else {
            if(a_qual[ia]>b_qual[ib]) {
                a_qual[ia] += 0.2*a_qual[ia];
                b_qual[ib] = 0;
            } else {
                b_qual[ib] += 0.2*b_qual[ib];
                a_qual[ia] = 0;
            }
        }
        a_qual[ia] = (a_qual[ia]<=255) ? a_qual[ia] : 255;
        b_qual[ib] = (b_qual[ib]<=255) ? b_qual[ib] : 255;
        ia++;
        ib++;
    }

quit :
    free(posa);
    free(posb);
}

int custom_overlap_constructor(void *data, const bam1_t *b, bam_pileup_cd *cd) {
    int ret;
    mplp_data *foo = (mplp_data*) data;
    khash_t(olap_hash) * ohash = ( khash_t(olap_hash) *)foo->ohash;
    khiter_t k = kh_get(olap_hash, ohash, bam_get_qname(b));
    bam1_t *a;
    // Skip unpaired reads
    if(!(b->core.flag & BAM_FPAIRED) || ((b->core.flag & 12) > 0)) return 0;
    if (k==kh_end(ohash)) {
        k = kh_put(olap_hash, ohash, bam_get_qname(b), &ret);
        kh_value(ohash, k) = b;
    } else {
        a = kh_value(ohash, k);
        cust_tweak_overlap_quality(a, b);
        kh_del(olap_hash, ohash, k);
    }
    return 0;
    
}

int custom_overlap_destructor(void *data, const bam1_t *b, bam_pileup_cd *cd) {
    mplp_data *foo = (mplp_data*) data;
    khash_t(olap_hash) * ohash = ( khash_t(olap_hash) *)foo->ohash;
    khiter_t k = kh_get(olap_hash, ohash, bam_get_qname(b));
    if (k!=kh_end(ohash)) kh_del(olap_hash, ohash, k);
    return 0;
}
