#include <inttypes.h>
#include <ctype.h>
#include <stdio.h>
#include <assert.h>
#include "htslib/khash.h"
#include "htslib/sam.h"

int getStrand(bam1_t *b);

/**********************************************
* Everything below here needs to track htslib!
**********************************************/
typedef struct {
    int k, x, y, end;
} cstate_t;

static cstate_t g_cstate_null = { -1, 0, 0, 0 };

typedef struct __linkbuf_t {
    bam1_t b;
    int32_t beg, end;
    cstate_t s;
    struct __linkbuf_t *next;
} lbnode_t;

typedef struct {
    int cnt, n, max;
    lbnode_t **buf;
} mempool_t;

// Dictionary of overlapping reads
KHASH_MAP_INIT_STR(olap_hash, lbnode_t *)
typedef khash_t(olap_hash) olap_hash_t;

struct __bam_plp_t {
    mempool_t *mp;
    lbnode_t *head, *tail, *dummy;
    int32_t tid, pos, max_tid, max_pos;
    int is_eof, max_plp, error, maxcnt;
    uint64_t id;
    bam_pileup1_t *plp;
    // for the "auto" interface only
    bam1_t *b;
    bam_plp_auto_f func;
    void *data;
    olap_hash_t *overlaps;
};

struct __bam_mplp_t {
    int n;
    uint64_t min, *pos;
    bam_plp_t *iter;
    int *n_plp;
    const bam_pileup1_t **plp;
};

static inline lbnode_t *mp_alloc(mempool_t *mp)
{
    ++mp->cnt;
    if (mp->n == 0) return (lbnode_t*)calloc(1, sizeof(lbnode_t));
    else return mp->buf[--mp->n];
}

/**********************************************
* Done reinventing the wheel
**********************************************/
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
            } else if(op == 1 || op == 4 || op == 5) { //I, S, H
                *(positions+offset) = -1;
                offset++;
            } else if(op == 2 || op == 3) { //D, N
                previous_position++;
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
    while(posa[ia]<0 && ia<na) ia++;
    while(posa[ib]<0 && ib<nb) ib++;
    if(ia==na || ib==nb) goto quit;

    //Go to the first overlapping position
    if(posa[ia]<posb[ib]) {
        while(posa[ia]<posb[ib] && ia<na) ia++;
    } else {
        while(posb[ib]<posa[ia] && ib<nb) ib++;
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

//This is a modified version of overlap_push()
static void cust_overlap_push(bam_plp_t iter, lbnode_t *node)
{
    if ( !iter->overlaps ) return;

    // mapped mates only
    if ( node->b.core.flag&BAM_FMUNMAP ) return;

    //Same chromosomes only
    if((&node->b)->core.tid != (&node->b)->core.mtid) return;

    //Overlap impossible
    if((&node->b)->core.pos < (&node->b)->core.mpos) {
        if(bam_endpos(&node->b) < (&node->b)->core.mpos) return;
    }

    khiter_t kitr = kh_get(olap_hash, iter->overlaps, bam_get_qname(&node->b));
    if ( kitr==kh_end(iter->overlaps) )
    {
        if((&node->b)->core.pos > (&node->b)->core.mpos) return; //The mate should be here, but isn't
        int ret;
        kitr = kh_put(olap_hash, iter->overlaps, bam_get_qname(&node->b), &ret);
        kh_value(iter->overlaps, kitr) = node;
    }
    else
    {
        lbnode_t *a = kh_value(iter->overlaps, kitr);
        cust_tweak_overlap_quality(&a->b, &node->b);
        kh_del(olap_hash, iter->overlaps, kitr);
        assert(a->end-1 == a->s.end);
        a->end = a->b.core.pos + bam_cigar2rlen(a->b.core.n_cigar, bam_get_cigar(&a->b));
        a->s.end = a->end - 1;
    }
}

//Essentially overlap_remove()
static void cust_overlap_remove(bam_plp_t iter, const bam1_t *b) {
    if(!iter->overlaps) return;

    khiter_t kitr;
    if(b) {
        kitr = kh_get(olap_hash, iter->overlaps, bam_get_qname(b));
        if ( kitr!=kh_end(iter->overlaps) )
            kh_del(olap_hash, iter->overlaps, kitr);
    } else {
        // remove all
        for (kitr = kh_begin(iter->overlaps); kitr<kh_end(iter->overlaps); kitr++)
            if ( kh_exist(iter->overlaps, kitr) ) kh_del(olap_hash, iter->overlaps, kitr);
    }
}


//This is essentially just bam_plp_push()
int cust_plp_push(bam_plp_t iter, const bam1_t *b)
{
    if (iter->error) return -1;
    if (b) {
        if (b->core.tid < 0) { cust_overlap_remove(iter, b); return 0; }
        // Skip only unmapped reads here, any additional filtering must be done in iter->func
        if (b->core.flag & BAM_FUNMAP) { cust_overlap_remove(iter, b); return 0; }
        if (iter->tid == b->core.tid && iter->pos == b->core.pos && iter->mp->cnt > iter->maxcnt)
        {
            cust_overlap_remove(iter, b);
            return 0;
        }
        bam_copy1(&iter->tail->b, b);
        cust_overlap_push(iter, iter->tail);
        iter->tail->beg = b->core.pos;
        iter->tail->end = b->core.pos + bam_cigar2rlen(b->core.n_cigar, bam_get_cigar(b));
        iter->tail->s = g_cstate_null; iter->tail->s.end = iter->tail->end - 1; // initialize cstate_t
        if (b->core.tid < iter->max_tid) {
            fprintf(stderr, "[cust_plp_push] the input is not sorted (chromosomes out of order)\n");
            iter->error = 1;
            return -1;
        }
        if ((b->core.tid == iter->max_tid) && (iter->tail->beg < iter->max_pos)) {
            fprintf(stderr, "[cust_plp_push] the input is not sorted (reads out of order)\n");
            iter->error = 1;
            return -1;
        }
        iter->max_tid = b->core.tid; iter->max_pos = iter->tail->beg;
        if (iter->tail->end > iter->pos || iter->tail->b.core.tid > iter->tid) {
            iter->tail->next = mp_alloc(iter->mp);
            iter->tail = iter->tail->next;
        }
    } else iter->is_eof = 1;
    return 0;
}

//This is essentially just bam_plp_auto
const bam_pileup1_t *cust_plp_auto(bam_plp_t iter, int *_tid, int *_pos, int *_n_plp)
{
    const bam_pileup1_t *plp;
    if (iter->func == 0 || iter->error) { *_n_plp = -1; return 0; }
    if ((plp = bam_plp_next(iter, _tid, _pos, _n_plp)) != 0) return plp;
    else { // no pileup line can be obtained; read alignments
        *_n_plp = 0;
        if (iter->is_eof) return 0;
        int ret;
        while ( (ret=iter->func(iter->data, iter->b)) >= 0) {
            if (cust_plp_push(iter, iter->b) < 0) {
                *_n_plp = -1;
                return 0;
            }
            if ((plp = bam_plp_next(iter, _tid, _pos, _n_plp)) != 0) return plp;
            // otherwise no pileup line can be returned; read the next alignment.
        }
        if ( ret < -1 ) { iter->error = ret; *_n_plp = -1; return 0; }
        bam_plp_push(iter, 0);
        if ((plp = bam_plp_next(iter, _tid, _pos, _n_plp)) != 0) return plp;
        return 0;
    }
}


//This is essentially just bam_mplp_auto from htslib
int cust_mplp_auto(bam_mplp_t iter, int *_tid, int *_pos, int *n_plp, const bam_pileup1_t **plp)
{
    int i, ret = 0;
    uint64_t new_min = (uint64_t)-1;
    for (i = 0; i < iter->n; ++i) {
        if (iter->pos[i] == iter->min) {
            int tid, pos;
            iter->plp[i] = cust_plp_auto(iter->iter[i], &tid, &pos, &iter->n_plp[i]);
            if ( iter->iter[i]->error ) return -1;
            iter->pos[i] = iter->plp[i] ? (uint64_t)tid<<32 | pos : 0;
        }
        if (iter->plp[i] && iter->pos[i] < new_min) new_min = iter->pos[i];
    }
    iter->min = new_min;
    if (new_min == (uint64_t)-1) return 0;
    *_tid = new_min>>32; *_pos = (uint32_t)new_min;
    for (i = 0; i < iter->n; ++i) {
        if (iter->pos[i] == iter->min) { // FIXME: valgrind reports "uninitialised value(s) at this line"
            n_plp[i] = iter->n_plp[i], plp[i] = iter->plp[i];
            ++ret;
        } else n_plp[i] = 0, plp[i] = 0;
    }
    return ret;
}
