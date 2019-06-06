#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include "MethylDackel.h"

void print_version(void);

extern inline double logit(double p) { 
    return(log(p) - log(1 - p)); 
}

//N.B., a tid of -1 means that the lastCall was written
struct lastCall{
    int32_t tid, pos;
    uint32_t nmethyl, nunmethyl;
};

const char *quux[2] = {"foo", "bar"};
const char *TriNucleotideContexts[25] = {"CAA", "CAC", "CAG", "CAT", "CAN", \
                                         "CCA", "CCC", "CCG", "CCT", "CCN", \
                                         "CGA", "CGC", "CGG", "CGT", "CGN", \
                                         "CTA", "CTC", "CTG", "CTT", "CTN", \
                                         "CNA", "CNC", "CNG", "CNT", "CNN"};

void writeCall(kstring_t *ks, Config *config, char *chrom, int32_t pos, int32_t width, uint32_t nmethyl, uint32_t nunmethyl, char base, char *context, const char *tnc) { 
    char str[10000]; // I don't really like hardcoding it, but given the probability that it ever won't suffice...
    char strand = (base=='C' || base=='c') ? 'F' : 'R';
    if(nmethyl+nunmethyl < config->minDepth && !config->cytosine_report) return;

    if (!config->fraction && !config->logit && !config->counts && !config->methylKit && !config->cytosine_report) {
        snprintf(str, 10000, \
            "%s\t%i\t%i\t%i\t%" PRIu32 "\t%" PRIu32 "\n", \
            chrom, \
            pos, \
            pos+width, \
            (int) (100.0*((double) nmethyl)/(nmethyl+nunmethyl)),\
            nmethyl, \
            nunmethyl);
    } else if(config->fraction) {
        snprintf(str, 10000, \
            "%s\t%i\t%i\t%f\n", \
            chrom, \
            pos, \
            pos+width, \
            ((double) nmethyl)/(nmethyl+nunmethyl));
    } else if(config->counts) {
        snprintf(str, 10000, \
            "%s\t%i\t%i\t%i\n", \
            chrom, \
            pos, \
            pos+width, \
            nmethyl+nunmethyl);
    } else if(config->logit) {
        snprintf(str, 10000, \
            "%s\t%i\t%i\t%f\n", \
            chrom, \
            pos, \
            pos+width, \
            logit(((double) nmethyl)/(nmethyl+nunmethyl)));
    } else if(config->methylKit) {
        snprintf(str, 10000, \
            "%s.%i\t%s\t%i\t%c\t%i\t%6.2f\t%6.2f\n", \
            chrom, \
            pos+1, \
            chrom, \
            pos+1, \
            strand, \
            nmethyl+nunmethyl, \
            100.0 * ((double) nmethyl)/(nmethyl+nunmethyl), \
            100.0 * ((double) nunmethyl)/(nmethyl+nunmethyl));
    } else if(config->cytosine_report) {
        strand = (base=='C' || base=='c') ? '+' : '-';
        snprintf(str, 10000, \
            "%s\t%i\t%c\t%"PRIu32"\t%"PRIu32"\tC%s\t%s\n", \
            chrom, \
            pos+1, \
            strand, \
            nmethyl, \
            nunmethyl, \
            context, \
            tnc);
    }

    kputs(str, ks);
}

char revcomp(char b) {
    switch(b) {
        case 'A':
        case 'a':
            return 'T';
        case 'C':
        case 'c':
            return 'G';
        case 'G':
        case 'g':
            return 'C';
        case 'T':
        case 't':
            return 'A';

            return 'N';
    }
}

int getTriNucContext(char *seq, uint32_t offset, int seqlen, int direction) {
    int rv = 0;
    char base;

    //Last base: column
    if((direction > 0 && offset+2 >= seqlen) || (direction < 0 && offset <= 1)) rv = 4;
    else {
        base = *(seq+offset+2*direction);
        if(direction < 0) base = revcomp(base);
        switch(base) {
            case 'A':
            case 'a':
                rv = 0;
                break;
            case 'C':
            case 'c':
                rv = 1;
                break;
            case 'G':
            case 'g':
                rv = 2;
                break;
            case 'T':
            case 't':
                rv = 3;
                break;
            default:
                rv = 4;
                break;
        }
    }

    //Middle
    if((direction > 0 && offset+1 >= seqlen) || (direction < 0 && offset == 0)) rv += 20;
    else {
        base = *(seq+offset+direction);
        if(direction < 0) base = revcomp(base);
        switch(base) {
            case 'A':
            case 'a':
                rv += 0;
                break;
            case 'C':
            case 'c':
                rv += 5;
                break;
            case 'G':
            case 'g':
                rv += 10;
                break;
            case 'T':
            case 't':
                rv += 15;
                break;
            default:
                rv += 20;
                break;
        }
    }
    return rv;
}

void writeBlank(kstring_t **ks, Config *config, char *chrom, int32_t pos, uint32_t localPos2, uint32_t *lastPos, char *seq, int seqlen) {
    int triNucContext = 0;
    int direction = 0;
    char context[3] = "HG";
    if(pos == -1) return;
    for(;*lastPos < pos; (*lastPos)++) {
        if((direction = isCpG(seq, *lastPos-localPos2, seqlen)) != 0) {
            if(!config->keepCpG) continue;
            triNucContext = getTriNucContext(seq, *lastPos - localPos2, seqlen, direction);
            context[0] = 'G'; context[1] = 0;
        } else if((direction = isCHG(seq, *lastPos-localPos2, seqlen)) != 0) {
            if(!config->keepCHG) continue;
            triNucContext = getTriNucContext(seq, *lastPos - localPos2, seqlen, direction);
            context[0] = 'H'; context[1] = 'G';
        } else if((direction = isCHH(seq, *lastPos-localPos2, seqlen)) != 0) {
            if(!config->keepCHH) continue;
            triNucContext = getTriNucContext(seq, *lastPos - localPos2, seqlen, direction);
            context[0] = 'H'; context[1] = 'H';
        } else {
            continue;
        }
        writeCall(ks[0], config, chrom, *lastPos, 1, 0, 0, (direction>0)?'C':'G', context, TriNucleotideContexts[triNucContext]);
    }
}

void processLast(kstring_t *ks, Config *config, struct lastCall *last, bam_hdr_t *hdr, int32_t tid, int32_t pos, int width, uint32_t nmethyl, uint32_t nunmethyl, char base) {
    if(last->tid == tid && last->pos == pos) {
        nmethyl += last->nmethyl;
        nunmethyl += last->nunmethyl;
        writeCall(ks, config, hdr->target_name[tid], pos, width, nmethyl, nunmethyl, base, NULL, NULL);
        last->tid = -1;
    } else {
        if(last->tid != -1) {
            writeCall(ks, config, hdr->target_name[last->tid], last->pos, width, last->nmethyl, last->nunmethyl, base, NULL, NULL);
        }
        last->tid = tid;
        last->pos = pos;
        last->nmethyl = nmethyl;
        last->nunmethyl = nunmethyl;
    }
}

// The opposite strand of a C should be a G. Ns are ignored
int isVariant(Config *config, const bam_pileup1_t *plp, uint32_t *coverage, int strand) {
    uint8_t base = bam_seqi(bam_get_seq(plp->b), plp->qpos);

    //Is the phred score even high enough?
    if(bam_get_qual(plp->b)[plp->qpos] < config->minPhred) return 0;

    *coverage += 1;
    if(strand & 1) { //OT or CTOT
        if(base != 4 && base != 15) return 1;
        else return 0;
    } else { // OB or CTOB
        if(base != 2 && base != 15) return 1;
        else return 0;
    }
}

void *extractCalls(void *foo) {
    //fprintf(stderr, "in extractCalls\n");
    Config *config = (Config*) foo;
    bam_hdr_t *hdr;
    //int counter_dbg = 0;
    bam_mplp_t iter;
    int ret, tid, pos, i, seqlen, type, rv, o = 0;
    int32_t bedIdx = 0;
    int n_plp; //This will need to be modified for multiple input files
    int strand, direction;
    uint32_t nmethyl = 0, nunmethyl = 0, nOff = 0, nVariant = 0;
    uint32_t localBin = 0, localPos = 0, localEnd = 0, localTid = 0, localPos2 = 0, lastPos = 0;
    uint64_t nVariantPositions = 0;
    const bam_pileup1_t **plp = NULL;
    char *seq = NULL, base = 'A';
    char context[3] = "HG";
    int tnc;
    mplp_data *data = NULL;
    struct lastCall *lastCpG = NULL;
    struct lastCall *lastCHG = NULL;
    kstring_t **os = NULL, *os_CpG = NULL, *os_CHG = NULL, *os_CHH = NULL;
    faidx_t *fai;
    hts_idx_t *bai;
    htsFile *fp;
    os = calloc(3, sizeof(kstring_t*));
    os_CpG = calloc(1, sizeof(kstring_t));
    os_CHG = calloc(1, sizeof(kstring_t));
    os_CHH = calloc(1, sizeof(kstring_t));
    if(!os_CpG || !os_CHG || !os_CHH || !os) {
        fprintf(stderr, "Couldn't allocate space for the kstring_t structures in extractCalls()!\n");
        return NULL;
    }
    os[0] = os_CpG;
    os[1] = os_CHG;
    os[2] = os_CHH;

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
    if(config->merge) {
        if(config->keepCpG) {
            lastCpG = calloc(1, sizeof(struct lastCall));
            assert(lastCpG);
            lastCpG->tid = -1;
        }
        if(config->keepCHG) {
            lastCHG = calloc(1, sizeof(struct lastCall));
            assert(lastCHG);
            lastCHG->tid = -1;
        }
    }

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
        //fprintf(stderr, "in extractCalls big loop\n");
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
        fprintf(stderr, "chrom name: %s\n", data->hdr->target_name[localTid]);
        //fprintf(stderr, "in extractCalls big loop (2)\n");
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
        //fprintf(stderr, "in extractCalls big loop (3)\n");
        localPos2 = 0;
        if(localPos > 1) {
            localPos2 = localPos - 2;
        }
        lastPos = localPos;


        //Break out of the loop if finished
        if(localTid >= hdr->n_targets) break; //Finish looping
        if(globalEnd && localPos >= globalEnd) break;
        data->iter = sam_itr_queryi(bai, localTid, localPos, localEnd);

        seq = faidx_fetch_seq(fai, hdr->target_name[localTid], localPos2, localEnd+10, &seqlen);
        if(seqlen < 0) {
            fprintf(stderr, "faidx_fetch_seq returned %i while trying to fetch the sequence for tid %s:%"PRIu32"-%"PRIu32"!\n",\
                seqlen, hdr->target_name[localTid], localPos2, localEnd);
            fprintf(stderr, "Note that the output will be truncated!\n");
            return NULL;
        }

        //Start the pileup
        //fprintf(stderr, "starting pileup\n");
        iter = bam_mplp_init(1, filter_func, (void **) &data);
        //fprintf(stderr, "starting pileup (2)\n");
        bam_mplp_init_overlaps(iter);
        //fprintf(stderr, "starting pileup (2b)\n");
        bam_mplp_set_maxcnt(iter, config->maxDepth);
        //fprintf(stderr, "starting pileup (2c)\n");
        while((ret = cust_mplp_auto(iter, &tid, &pos, &n_plp, plp)) > 0) {
            //fprintf(stderr, "looping %d\n", counter_dbg++);
            if(pos < localPos || pos >= localEnd) continue; // out of the region requested

            if(config->bed) { //Handle -l
                while((o = posOverlapsBED(tid, pos, config->bed, bedIdx)) == -1) bedIdx++;
                if(o == 0) continue; //Wrong strand
            }

            if((direction = isCpG(seq, pos-localPos2, seqlen))) {
                if(!config->keepCpG) continue;
                type = 0;
            } else if((direction = isCHG(seq, pos-localPos2, seqlen))) {
                if(!config->keepCHG) continue;
                type = 1;
            } else if((direction = isCHH(seq, pos-localPos2, seqlen))) {
                if(!config->keepCHH) continue;
                type = 2;
            } else {
                continue;
            }

            nmethyl = nunmethyl = nVariant = nOff = 0;
            base = *(seq+pos-localPos2);
            for(i=0; i<n_plp; i++) {
                if(plp[0][i].is_del) continue;
                if(plp[0][i].is_refskip) continue;
                if(config->bed) if(!readStrandOverlapsBED(plp[0][i].b, config->bed->region[bedIdx])) continue;
                strand = getStrand((plp[0]+i)->b);
                if(strand & 1) {
                    if(base != 'C' && base != 'c') {
                        nVariant += isVariant(config, plp[0]+i, &nOff, strand);
                        continue;
                    }
                } else {
                    if(base != 'G' && base != 'g') {
                        nVariant += isVariant(config, plp[0]+i, &nOff, strand);
                        continue;
                    }
                }
                rv = updateMetrics(config, plp[0]+i);
                if(rv > 0) nmethyl++;
                else if(rv<0) nunmethyl++;
            }

            // Skip likely variant positions
            if(config->minOppositeDepth > 0 && \
               nOff >= config->minOppositeDepth && \
               ((double) nVariant) / ((double) nOff) >= config->maxVariantFrac) {
                nVariantPositions++;
                //If we're merging context, then skip an entire merged site
                if(config->merge) {
                    if(type == 0 && lastCpG->tid == tid && lastCpG->pos == pos - 1 && (base == 'G' || base == 'g')) {
                        lastCpG->nmethyl = 0;
                        lastCpG->nunmethyl = 0;
                    } else if(type == 1 && lastCHG->tid == tid && lastCHG->pos == pos - 2 && (base == 'G' || base == 'g')) {
                        lastCHG->nmethyl = 0;
                        lastCHG->nunmethyl = 0;
                    }
                }
                continue;
            }

            if(nmethyl+nunmethyl==0 && config->cytosine_report == 0) continue;
            if(!config->merge || type==2) { //Also, cytosine report
                if(config->cytosine_report) {
                    writeBlank(os, config, hdr->target_name[tid], pos, localPos2, &lastPos, seq, seqlen);

                    //Set the C-context
                    if(type == 0) {
                        context[0] = 'G'; context[1] = 0;
                    } else if(type == 1) {
                        context[0] = 'H'; context[1] = 'G';
                    } else {
                        context[0] = 'H'; context[1] = 'H';
                    }

                    //Set the trinucleotide context
                    tnc = getTriNucContext(seq, pos - localPos2, seqlen, direction);

                    writeCall(os[0], config, hdr->target_name[tid], pos, 1, nmethyl, nunmethyl, base, context, TriNucleotideContexts[tnc]);
                } else {
                    writeCall(os[type], config, hdr->target_name[tid], pos, 1, nmethyl, nunmethyl, base, NULL, NULL);
                }
            } else {
                //Merge into per-CpG/CHG metrics
                if(type==0) {
                    if(base=='G' || base=='g') pos--;
                    processLast(os_CpG, config, lastCpG, hdr, tid, pos, 2, nmethyl, nunmethyl, base);
                } else {
                    if(base=='G' || base=='g') pos-=2;
                    processLast(os_CHG, config, lastCHG, hdr, tid, pos, 3, nmethyl, nunmethyl, base);
                }
            }
            lastPos = pos+1;
        }
        //fprintf(stderr, "starting pileup (3)\n");
        bam_mplp_destroy(iter);

        //Don't forget the last CpG/CHG
        nmethyl = 0;
        nunmethyl = 0;
        if(config->merge) {
            if(config->keepCpG && lastCpG->tid != -1) {
                processLast(os_CpG, config, lastCpG, hdr, tid, pos, 2, nmethyl, nunmethyl, base);
                lastCpG->tid = -1;
            }
            if(config->keepCHG && lastCHG->tid != -1) {
                processLast(os_CHG, config, lastCHG, hdr, tid, pos, 3, nmethyl, nunmethyl, base);
                lastCHG->tid = -1;
            }
        } else if(config->cytosine_report) {
            writeBlank(os, config, hdr->target_name[tid], localEnd, localPos2, &lastPos, seq, seqlen);
        }
        hts_itr_destroy(data->iter);
        free(seq);

        while(1) {
            pthread_mutex_lock(&outputMutex);
            if(outputBin != localBin) {
                pthread_mutex_unlock(&outputMutex);
                continue;
            }
            if(config->keepCpG && os_CpG->l) {
                fputs(os_CpG->s, config->output_fp[0]);
                os_CpG->l = 0;
            }
            if(config->keepCHG && os_CHG->l) {
                fputs(os_CHG->s, config->output_fp[1]);
                os_CHG->l = 0;
            }
            if(config->keepCHH && os_CHH->l) {
                fputs(os_CHH->s, config->output_fp[2]);
                os_CHH->l = 0;
            }
            outputBin++;
            pthread_mutex_unlock(&outputMutex);
            break;
        }
    }

    free(os_CpG->s); free(os_CpG);
    free(os_CHG->s); free(os_CHG);
    free(os_CHH->s); free(os_CHH);
    free(os);
    bam_hdr_destroy(hdr);
    fai_destroy(fai);
    hts_close(fp);
    hts_idx_destroy(bai);
    free(data);
    free(plp);
    if(config->merge) {
        if(config->keepCpG) free(lastCpG);
        if(config->keepCHG) free(lastCHG);
    }

    if(nVariantPositions > 0) {
        pthread_mutex_lock(&outputMutex);
        globalnVariantPositions += nVariantPositions;
        pthread_mutex_unlock(&outputMutex);
    }
    //fprintf(stderr, "done with extractCalls\n");
    return NULL;
}

void printHeader(FILE *of, char *context, char *opref, Config config) {
    fprintf(of, "track type=\"bedGraph\" description=\"%s %s", opref, context);
    if(config.merge) fprintf(of, " merged");
    if(config.fraction) fprintf(of, " methylation fractions\"\n");
    else if(config.counts) fprintf(of, " methylation counts\"\n");
    else if(config.logit) fprintf(of, " logit transformed methylation fractions\"\n");
    else fprintf(of, " methylation levels\"\n");
}

void extract_usage() {
    fprintf(stderr, "\nUsage: MethylDackel extract [OPTIONS] <ref.fa> <sorted_alignments.bam>\n");
    fprintf(stderr,
"\n"
"Options:\n"
" -q INT           Minimum MAPQ threshold to include an alignment (default 10)\n"
" -p INT           Minimum Phred threshold to include a base (default 5). This\n"
"                  must be >0.\n"
" -D INT           Maximum per-base depth (default 2000)\n"
" -d INT           Minimum per-base depth for reporting output. If you use\n"
"                  --mergeContext, this then applies to the merged CpG/CHG.\n"
"                  (default 1)\n"
" -r STR           Region string in which to extract methylation\n"
" -l FILE          A BED file listing regions for inclusion.\n"
" --keepStrand     If a BED file is specified, then this option will cause the\n"
"                  strand column (column 6) to be utilized, if present. Thus, if\n"
"                  a region has a '+' in this column, then only metrics from the\n"
"                  top strand will be output. Note that the -r option can be used\n"
"                  to limit the regions of -l.\n"
" -M, --mappability FILE    A bigWig file containing mappability data for\n"
"                  filtering reads.\n"
" -t, --mappabilityThreshold FLOAT    If a bigWig file is provided, this sets the\n"
"                  threshold mappability value above which a base is considered\n"
"                  mappable (default 0.01).\n"
" -b, --minMappableBases INT    If a bigWig file is provided, this sets the\n"
"                  number of mappable bases needed for a read to be considered\n"
"                  mappable (default 15).\n"
" -@ nThreads      The number of threads to use, the default 1\n"
" --chunkSize INT  The size of the genome processed by a single thread at a time.\n"
"                  The default is 1000000 bases. This value MUST be at least 1.\n"
" --mergeContext   Merge per-Cytosine metrics from CpG and CHG contexts into\n"
"                  per-CPG or per-CHG metrics.\n"
" -o, --opref STR  Output filename prefix. CpG/CHG/CHH context metrics will be\n"
"                  output to STR_CpG.bedGraph and so on.\n"
" --keepDupes      By default, any alignment marked as a duplicate is ignored.\n"
"                  This option causes them to be incorporated.\n"
" --keepSingleton  By default, if only one read in a pair aligns (a singleton)\n"
"                  then it's ignored.\n"
" --keepDiscordant By default, paired-end alignments with the properly-paired bit\n"
"                  unset in the FLAG field are ignored. Note that the definition\n"
"                  of concordant and discordant is based on your aligner\n"
"                  settings.\n"
" -F, --ignoreFlags    By default, any alignment marked as secondary (bit 0x100),\n"
"                  failing QC (bit 0x200), a PCR/optical duplicate (0x400) or\n"
"                  supplemental (0x800) is ignored. This equates to a value of\n"
"                  0xF00 or 3840 in decimal. If you would like to change that,\n"
"                  you can specify a new value here.\n"
"                  ignored. Specifying this causes them to be included.\n"
" -R, --requireFlags   Require each alignment to have all bits in this value\n"
"                  present, or else the alignment is ignored. This is equivalent\n"
"                  to the -f option in samtools. The default is 0, which\n"
"                  includes all alignments.\n"
" --noCpG          Do not output CpG context methylation metrics\n"
" --CHG            Output CHG context methylation metrics\n"
" --CHH            Output CHH context methylation metrics\n"
" --fraction       Extract fractional methylation (only) at each position. This\n"
"                  produces a file with a .meth.bedGraph extension.\n"
" --counts         Extract base counts (only) at each position. This produces a\n"
"                  file with a .counts.bedGraph extension.\n"
" --logit          Extract logit(M/(M+U)) (only) at each position. This produces\n"
"                  a file with a .logit.bedGraph extension.\n"
" --minOppositeDepth   If you would like to exclude sites that likely contain\n"
"                  SNPs, then specifying a value greater than 0 here will\n"
"                  indicate the minimum depth required on the strand opposite of\n"
"                  a C to look for A/T/C bases. The fraction of these necessary\n"
"                  to exclude a position from methylation extraction is specified\n"
"                  by the --maxVariantFrac parameter. The default is 0, which\n"
"                  means that no positions will be excluded. Note that the -p and\n"
"                  -q apply here as well. Note further that if you use\n"
"                  --mergeContext that a merged site will be excluded if either\n"
"                  of its individual Cs would be excluded.\n"
" --maxVariantFrac The maximum fraction of A/T/C base calls on the strand\n"
"                  opposite of a C to allow before a position is declared a\n"
"                  variant (thereby being excluded from output). The default is\n"
"                  0.0. See also --minOppositeDepth.\n"
" --methylKit      Output in the format required by methylKit. Note that this is\n"
"                  incompatible with --mergeContext, --fraction and --counts.\n"
" --cytosine_report  A per-base exhaustive report comparable to that produced\n"
"                  with the same option in Bismark's methylation extractor. The\n"
"                  output is a tab-separated file with the following columns:\n"
"                  chromosome, position (1-based), strand, number of alignments\n"
"                  supporting methylation, number of alignments supporting\n"
"                  unmethylation, CG/CHG/CHH, trinucleotide context. This is not\n"
"                  compatible with --fraction, --counts, --methylKit, or\n"
"                  --mergeContext. The produces a single file with a\n"
"                  .cytosine_report.txt extension. Note that even bases with no\n"
"                  coverage will be included in the output.\n"
" --OT INT,INT,INT,INT Inclusion bounds for methylation calls from reads/pairs\n"
"                  originating from the original top strand. Suggested values can\n"
"                  be obtained from the MBias program. Each integer represents a\n"
"                  1-based position on a read. For example --OT A,B,C,D\n"
"                  translates to, \"Include calls at positions from A through B\n"
"                  on read #1 and C through D on read #2\". If a 0 is used a any\n"
"                  position then that is translated to mean start/end of the\n"
"                  alignment, as appropriate. For example, --OT 5,0,0,0 would\n"
"                  include all but the first 4 bases on read #1. Users are\n"
"                  strongly advised to consult a methylation bias plot, for\n"
"                  example by using the MBias program.\n"
" --OB INT,INT,INT,INT\n"
" --CTOT INT,INT,INT,INT\n"
" --CTOB INT,INT,INT,INT As with --OT, but for the original bottom, complementary\n"
"                  to the original top, and complementary to the original bottom\n"
"                  strands, respectively.\n"
" --nOT INT,INT,INT,INT Like --OT, but always exclude INT bases from a given end\n"
"                  from inclusion,regardless of the length of an alignment. This\n"
"                  is useful in cases where reads may have already been trimmed\n"
"                  to different lengths, but still none-the-less contain a\n"
"                  certain length bias at one or more ends.\n"
" --nOB INT,INT,INT,INT\n"
" --nCTOT INT,INT,INT,INT\n"
" --nCTOB INT,INT,INT,INT As with --nOT, but for the original bottom,\n"
"                  complementary to the original top, and complementary to the\n"
"                  original bottom strands, respectively.\n"
" --version        Print version and then quit.\n"
"\nNote that --fraction, --counts, and --logit are mutually exclusive!\n");
}

int extract_main(int argc, char *argv[]) {
    //fprintf(stderr, "in extract_main\n");
    char *opref = NULL, *oname, *p;
    int c, i, keepStrand = 0;
    Config config;
    bam_hdr_t *hdr = NULL;

    globalTid = globalPos = globalEnd = bin = globalnVariantPositions = 0;

    //Defaults
    config.BWName = NULL;
    config.BW_ptr = NULL;
    config.mappabilityCutoff = 0.01;
    config.minMappableBases = 15;
    config.keepCpG = 1; config.keepCHG = 0; config.keepCHH = 0;
    config.minMapq = 10; config.minPhred = 5; config.keepDupes = 0;
    config.keepSingleton = 0, config.keepDiscordant = 0;
    config.minDepth = 1;
    config.methylKit = 0;
    config.merge = 0;
    config.minOppositeDepth = 0;
    config.maxVariantFrac = 0.0;
    config.maxDepth = 2000;
    config.fp = NULL;
    config.bai = NULL;
    config.reg = NULL;
    config.bedName = NULL;
    config.bed = NULL;
    config.fraction = 0;
    config.counts = 0;
    config.logit = 0;
    config.ignoreFlags = 0xF00;
    config.requireFlags = 0;
    config.nThreads = 1;
    config.chunkSize = 1000000;
    config.cytosine_report = 0;
    for(i=0; i<16; i++) config.bounds[i] = 0;
    for(i=0; i<16; i++) config.absoluteBounds[i] = 0;

    static struct option lopts[] = {
        {"opref",        1, NULL, 'o'},
        {"fraction",     0, NULL, 'f'},
        {"counts",       0, NULL, 'c'},
        {"logit",        0, NULL, 'm'},
        {"minDepth",     1, NULL, 'd'},
        {"noCpG",        0, NULL,   1},
        {"CHG",          0, NULL,   2},
        {"CHH",          0, NULL,   3},
        {"keepDupes",    0, NULL,   4},
        {"keepSingleton",0, NULL,   5},
        {"keepDiscordant",0,NULL,   6},
        {"OT",           1, NULL,   7},
        {"OB",           1, NULL,   8},
        {"CTOT",         1, NULL,   9},
        {"CTOB",         1, NULL,  10},
        {"mergeContext", 0, NULL,  11},
        {"methylKit",    0, NULL,  12},
        {"nOT",          1, NULL,  13},
        {"nOB",          1, NULL,  14},
        {"nCTOT",        1, NULL,  15},
        {"nCTOB",        1, NULL,  16},
        {"minOppositeDepth", 1, NULL, 17},
        {"maxVariantFrac", 1, NULL, 18},
        {"chunkSize",    1, NULL,  19},
        {"keepStrand",   0, NULL,  20},
        {"cytosine_report", 0, NULL, 21},
        {"ignoreFlags",  1, NULL, 'F'},
        {"requireFlags", 1, NULL, 'R'},
        {"help",         0, NULL, 'h'},
        {"version",      0, NULL, 'v'},
        {"mappability",         1, NULL,  'M'},
        {"mappabilityThreshold",         1, NULL,  't'},
        {"minMappableBases",         1, NULL,  'b'},
        {0,              0, NULL,   0}
    };
    while((c = getopt_long(argc, argv, "hvq:p:r:l:o:D:f:c:m:d:F:R:@:M:t:b:", lopts,NULL)) >=0){
        switch(c) {
        case 'h':
            extract_usage();
            return 0;
        case 'v':
            print_version();
            return 0;
        case 'o':
            opref = strdup(optarg);
            break;
        case 'D':
            config.maxDepth = atoi(optarg);
            break;
        case 'd':
            config.minDepth = atoi(optarg);
            if(config.minDepth < 1) {
                fprintf(stderr, "Error, the minimum depth must be at least 1!\n");
                return 1;
            }
            break;
	case 'r':
	    config.reg = optarg;
	    break;
        case 'l':
            config.bedName = optarg;
            break;
        case 1:
            config.keepCpG = 0;
            break;
        case 2:
            config.keepCHG = 1;
            break;
        case 3:
            config.keepCHH = 1;
            break;
        case 4:
            config.keepDupes = 1;
            break;
        case 5:
            config.keepSingleton = 1;
            break;
        case 6:
            config.keepDiscordant = 1;
            break;
        case 7:
            parseBounds(optarg, config.bounds, 0);
            break;
        case 8:
            parseBounds(optarg, config.bounds, 1);
            break;
        case 9:
            parseBounds(optarg, config.bounds, 2);
            break;
        case 10:
            parseBounds(optarg, config.bounds, 3);
            break;
        case 11:
            config.merge = 1;
            break;
        case 12:
            config.methylKit = 1;
            break;
        case 13:
            parseBounds(optarg, config.absoluteBounds, 0);
            break;
        case 14:
            parseBounds(optarg, config.absoluteBounds, 1);
            break;
        case 15:
            parseBounds(optarg, config.absoluteBounds, 2);
            break;
        case 16:
            parseBounds(optarg, config.absoluteBounds, 3);
            break;
        case 17:
            config.minOppositeDepth = atoi(optarg);
            break;
        case 18:
            config.maxVariantFrac = atof(optarg);
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
            config.cytosine_report = 1;
            break;
        case 'M':
            config.BWName = optarg;
            break;
        case 't':
            config.mappabilityCutoff = atof(optarg);
            break;
        case 'b':
            config.minMappableBases = atoi(optarg);
            break;
        case 'F':
            config.ignoreFlags = atoi(optarg);
            break;
        case 'R':
            config.requireFlags = atoi(optarg);
            break;
        case 'q':
            config.minMapq = atoi(optarg);
            break;
        case 'p':
            config.minPhred = atoi(optarg);
            break;
        case 'm':
            config.logit = 1;
            break;
        case 'f':
            config.fraction = 1;
            break;
        case 'c':
            config.counts = 1;
            break;
        case '@':
            config.nThreads = atoi(optarg);
            break;
        case '?':
        default :
            fprintf(stderr, "Invalid option '%c'\n", c);
            extract_usage();
            return 1;
        }
    }

    if(argc == 1) {
        extract_usage();
        return 0;
    }
    if(argc-optind != 2) {
        fprintf(stderr, "You must supply a reference genome in fasta format and an input BAM file!!!\n");
        extract_usage();
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
    if(config.fraction+config.counts+config.logit+config.methylKit+config.cytosine_report > 1) {
        fprintf(stderr, "More than one of --fraction, --counts, --methylKit, --cytosine_report and --logit were specified. These are mutually exclusive.\n");
        extract_usage();
        return 1;
    }
    if(config.methylKit + config.merge == 2) {
        fprintf(stderr, "--mergeContext and --methylKit are mutually exclusive.\n");
        extract_usage();
        return 1;
    }
    if(config.cytosine_report + config.merge == 2) {
        fprintf(stderr, "--mergeContext and --cytosine_report are mutually exclusive.\n");
        extract_usage();
        return 1;
    }

    //Has more than one output format been requested?
    if(config.fraction + config.counts + config.logit > 1) {
        fprintf(stderr, "You may specify AT MOST one of -c/--counts, -f/--fraction, or -m/--logit.\n");
        return -6;
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
    if(config.BWName && (config.BW_ptr = bwOpen(config.BWName, NULL, "r")) == NULL) {
        fprintf(stderr, "Couldn't open %s for reading!\n", config.BWName);
        return -4;
    }
    if(config.BWName)
    {
        double lastVal = 0;
        config.bw_data = malloc(config.BW_ptr->cl->nKeys*sizeof(char*)); //init outer array
        fprintf(stderr, "loading mappability data from %s\n", config.BWName);
        //fprintf(stderr, "chrom, len\n");
        for(int i = 0; i<config.BW_ptr->cl->nKeys; i++)
        {
            int arrlen;
            //fprintf(stderr, "%s: %d\n", config.BW_ptr->cl->chrom[i], config.BW_ptr->cl->len[i]);
            arrlen = config.BW_ptr->cl->len[i]/8;
            if(config.BW_ptr->cl->len[i]%8 > 0)
            {
                arrlen++;
            }
            config.bw_data[i] = malloc(arrlen*sizeof(char)); //init inner array
            //fprintf(stderr, "getting values...\n");
            bwOverlappingIntervals_t *vals = bwGetValues(config.BW_ptr, config.BW_ptr->cl->chrom[i], 0, config.BW_ptr->cl->len[i], 1);
            //fprintf(stderr, "loaded values, saving...\n");
            for(int j = 0; j<config.BW_ptr->cl->len[i]; j++)
            {
                char offset;
                int index;
                char aboveCutoff;
                double val;
                index = j/8;
                offset = j%8;
                //fprintf(stderr, "index: %d, offset: %d\n", index, offset);
                if(offset == 0) //starting new byte
                {
                    config.bw_data[i][index] = 0; //init new byte
                }
                val = vals->value[j];
                if(isnan(val))
                {
                    val = lastVal; //convert NA to previous value
                }
                else
                {
                    //fprintf(stderr, "new value: %f\n", val);
                    lastVal = val; //set new previous value
                }
                aboveCutoff = (char)(val > config.mappabilityCutoff); //check if above cutoff
                //fprintf(stderr, "val: %f, cutoff: %f, aboveCutoff: %d\n", val, config.mappabilityCutoff, aboveCutoff);
                config.bw_data[i][index] = config.bw_data[i][index] | (aboveCutoff << offset); //set bit
                //fprintf(stderr, "byte: 0x%x\n", config.bw_data[i][index]);
            }
            bwDestroyOverlappingIntervals(vals);
            
        }
    }

    //Output files
    //fprintf(stderr, "setting up output file\n");
    config.output_fp = malloc(sizeof(FILE *) * 3);
    assert(config.output_fp);
    if(opref == NULL) {
        opref = strdup(argv[optind+1]);
        assert(opref);
        p = strrchr(opref, '.');
        if(p != NULL) *p = '\0';
        fprintf(stderr, "writing to prefix:'%s'\n", opref);
        //fprintf(stderr, "printed \"writing to prefix\" message\n");
    }
    //fprintf(stderr, "starting processing\n");
    if(config.fraction) { 
        oname = malloc(sizeof(char) * (strlen(opref)+19));
    } else if(config.counts) {
        oname = malloc(sizeof(char) * (strlen(opref)+21));
    } else if(config.logit) {
        oname = malloc(sizeof(char) * (strlen(opref)+20));
    } else if(config.methylKit) {
        oname = malloc(sizeof(char) * (strlen(opref)+15));
    } else if(config.cytosine_report) {
        oname = malloc(sizeof(char) * (strlen(opref)+21));
        sprintf(oname, "%s.cytosine_report.txt", opref);
        config.output_fp[0] = fopen(oname, "w");
        config.output_fp[1] = config.output_fp[0];
        config.output_fp[2] = config.output_fp[0];
    } else { 
        oname = malloc(sizeof(char) * (strlen(opref)+14));
    }
    //fprintf(stderr, "read some config stuff\n");
    assert(oname);
    if(config.keepCpG && !config.cytosine_report) {
        if(config.fraction) { 
            sprintf(oname, "%s_CpG.meth.bedGraph", opref);
        } else if(config.counts) {
            sprintf(oname, "%s_CpG.counts.bedGraph", opref);
        } else if(config.logit) {
            sprintf(oname, "%s_CpG.logit.bedGraph", opref);
        } else if(config.methylKit) {
            sprintf(oname, "%s_CpG.methylKit", opref);
        } else { 
            sprintf(oname, "%s_CpG.bedGraph", opref);
        }
        config.output_fp[0] = fopen(oname, "w");
        if(config.output_fp[0] == NULL) {
            fprintf(stderr, "Couldn't open the output CpG metrics file for writing! Insufficient permissions?\n");
            return -3;
        }
        if(config.methylKit) {
            fprintf(config.output_fp[0], "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n");
        } else {
            printHeader(config.output_fp[0], "CpG", opref, config);
        }
    }
    //fprintf(stderr, "more config stuff\n");
    if(config.keepCHG && !config.cytosine_report) {
        if(config.fraction) { 
            sprintf(oname, "%s_CHG.meth.bedGraph", opref);
        } else if(config.counts) {
            sprintf(oname, "%s_CHG.counts.bedGraph", opref);
        } else if(config.logit) {
            sprintf(oname, "%s_CHG.logit.bedGraph", opref);
        } else if(config.methylKit) {
            sprintf(oname, "%s_CHG.methylKit", opref);
        } else { 
            sprintf(oname, "%s_CHG.bedGraph", opref);
        }
        config.output_fp[1] = fopen(oname, "w");
        if(config.output_fp[1] == NULL) {
            fprintf(stderr, "Couldn't open the output CHG metrics file for writing! Insufficient permissions?\n");
            return -3;
        }
        if(config.methylKit) {
            fprintf(config.output_fp[1], "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n");
        } else {
            printHeader(config.output_fp[1], "CHG", opref, config);
        }
    }
    //fprintf(stderr, "more config stuff (2)\n");
    if(config.keepCHH && !config.cytosine_report) {
        if(config.fraction) { 
            sprintf(oname, "%s_CHH.meth.bedGraph", opref);
        } else if(config.counts) {
            sprintf(oname, "%s_CHH.counts.bedGraph", opref);
        } else if(config.logit) {
            sprintf(oname, "%s_CHH.logit.bedGraph", opref);
        } else if(config.methylKit) {
            sprintf(oname, "%s_CHH.methylKit", opref);
        } else { 
            sprintf(oname, "%s_CHH.bedGraph", opref);
        }
        config.output_fp[2] = fopen(oname, "w");
        if(config.output_fp[2] == NULL) {
            fprintf(stderr, "Couldn't open the output CHH metrics file for writing! Insufficient permissions?\n");
            return -3;
        }
        if(config.methylKit) {
            fprintf(config.output_fp[2], "chrBase\tchr\tbase\tstrand\tcoverage\tfreqC\tfreqT\n");
        } else {
            printHeader(config.output_fp[2], "CHH", opref, config);
        }
    }
    //fprintf(stderr, "about to parse BED region (or not)\n");
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
    //fprintf(stderr, "starting pileup\n");
    //Run the pileup
    pthread_mutex_init(&positionMutex, NULL);
    pthread_mutex_init(&bwMutex, NULL);
    //fprintf(stderr, "starting threads\n");
    pthread_t *threads = calloc(config.nThreads, sizeof(pthread_t));
    //fprintf(stderr, "starting threads2\n");
    for(i=0; i < config.nThreads; i++) pthread_create(threads+i, NULL, &extractCalls, &config);
    //fprintf(stderr, "starting threads3\n");
    for(i=0; i < config.nThreads; i++) pthread_join(threads[i], NULL);
    //fprintf(stderr, "done with threads\n");
    free(threads);
    pthread_mutex_destroy(&bwMutex);
    pthread_mutex_destroy(&positionMutex);

    //If we've filtered out variant sites
    if(globalnVariantPositions) printf("%"PRIu64" positions were excluded due to likely being variants.\n", globalnVariantPositions);
    //Close things up
    hts_close(config.fp);
    if(config.cytosine_report) fclose(config.output_fp[0]);
    if(config.keepCpG && !config.cytosine_report) fclose(config.output_fp[0]);
    if(config.keepCHG && !config.cytosine_report) fclose(config.output_fp[1]);
    if(config.keepCHH && !config.cytosine_report) fclose(config.output_fp[2]);
    hts_idx_destroy(config.bai);
    free(opref);
    if(config.bed) destroyBED(config.bed);
    free(oname);
    free(config.output_fp);
    //fprintf(stderr, "there should be %d chromosomes\n", config.BW_ptr->cl->nKeys); 
    //fprintf(stderr, "config.bw_data[0] is \"%s\"\n", config.bw_data[0]);
    for(int i = 0; i<config.BW_ptr->cl->nKeys; i++)
    {
	//fprintf(stderr, "config.bw_data[%d] is \"%s\"\n", i, config.bw_data[i]);
        free(config.bw_data[i]);
    }
    free(config.bw_data);
    bwClose(config.BW_ptr);
    return 0;
}
