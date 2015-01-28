#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <assert.h>
#include <errno.h>
#include <limits.h>
#include "PileOMeth.h"

void extractCalls(Config *config) {
    bam_hdr_t *hdr = sam_hdr_read(config->fp);
    bam_mplp_t iter;
    int ret, tid, pos, i, seqlen, type, rv, o = 0;
    int beg0 = 0, end0 = 1u<<29;
    int n_plp; //This will need to be modified for multiple input files
    int ctid = -1; //The tid of the contig whose sequence is stored in "seq"
    int idxBED = 0, strand;
    uint32_t nmethyl, nunmethyl;
    const bam_pileup1_t **plp = NULL;
    char *seq = NULL, base;
    mplp_data *data = NULL;

    data = calloc(1,sizeof(mplp_data));
    if(data == NULL) {
        fprintf(stderr, "Couldn't allocate space for the data structure in extractCalls()!\n");
        return;
    }
    data->config = config;
    data->hdr = hdr;
    if (config->reg) {
        if((data->iter = sam_itr_querys(config->bai, hdr, config->reg)) == 0){
            fprintf(stderr, "failed to parse regions %s", config->reg);
            return;
        }
    }
    if(config->bedName) {
        config->bed = parseBED(config->bedName, hdr);
        if(config->bed == NULL) return;
    }

    plp = calloc(1, sizeof(bam_pileup1_t *)); //This will have to be modified for multiple input files
    if(plp == NULL) {
        fprintf(stderr, "Couldn't allocate space for the plp structure in extractCalls()!\n");
        return;
    }

    //Start the pileup
    iter = bam_mplp_init(1, filter_func, (void **) &data);
    bam_mplp_init_overlaps(iter);
    bam_mplp_set_maxcnt(iter, config->maxDepth);
    while((ret = cust_mplp_auto(iter, &tid, &pos, &n_plp, plp)) > 0) {
        //Do we need to process this position?
	if (config->reg){
	    beg0 = data->iter->beg, end0 = data->iter->end;
	    if ((pos < beg0 || pos >= end0)) continue; // out of the region requested
	}
        if(tid != ctid) {
            if(seq != NULL) free(seq);
            seq = faidx_fetch_seq(config->fai, hdr->target_name[tid], 0, faidx_seq_len(config->fai, hdr->target_name[tid]), &seqlen);
            if(seqlen < 0) {
                fprintf(stderr, "faidx_fetch_seq returned %i while trying to fetch the sequence for tid %i (%s)!\n",\
                    seqlen, tid, hdr->target_name[tid]);
                fprintf(stderr, "Note that the output will be truncated!\n");
                return;
            }
            ctid = tid;
        }

        if(config->bed) { //Handle -l
            while((o = posOverlapsBED(tid, pos, config->bed, idxBED)) == -1) idxBED++;
            if(o == 0) continue; //Wrong strand
        }

        if(isCpG(seq, pos, seqlen)) {
            if(!config->keepCpG) continue;
            type = 0;
        } else if(isCHG(seq, pos, seqlen)) {
            if(!config->keepCHG) continue;
            type = 1;
        } else if(isCHH(seq, pos, seqlen)) {
            if(!config->keepCHH) continue;
            type = 2;
        } else {
            continue;
        }

        nmethyl = nunmethyl = 0;
        base = *(seq+pos);
        for(i=0; i<n_plp; i++) {
            if(plp[0][i].is_del) continue;
            if(plp[0][i].is_refskip) continue;
            if(config->bed) if(!readStrandOverlapsBED(plp[0][i].b, config->bed->region[idxBED])) continue;
            strand = getStrand((plp[0]+i)->b);
            if(strand & 1) {
                if(base != 'C' && base != 'c') continue;
            } else {
                if(base != 'G' && base != 'g') continue;
            }
            rv = updateMetrics(config, plp[0]+i);
            if(rv > 0) nmethyl++;
            else if(rv<0) nunmethyl++;
        }

        if(nmethyl+nunmethyl) fprintf(config->output_fp[type], "%s\t%i\t%i\t%i\t%" PRIu32 "\t%" PRIu32 "\n", \
            hdr->target_name[tid], pos, pos+1, (int) (1000.0 * ((double) nmethyl)/(nmethyl+nunmethyl)), nmethyl, nunmethyl);
    }

    bam_hdr_destroy(hdr);
    if(data->iter) hts_itr_destroy(data->iter);
    bam_mplp_destroy(iter);
    free(data);
    free(plp);
    if(seq != NULL) free(seq);
}

void parseBounds(char *s2, int *vals, int mult) {
    char *p, *s = strdup(s2), *end;
    int i, v;

    p = strtok(s, ",");
    v = strtol(p, &end, 10);
    if((errno == ERANGE && (v == LONG_MAX || v == LONG_MIN)) || (errno != 0 && v == 0) || end == p) v = -1;
    if(v>=0) vals[4*mult] = v;
    else {
        fprintf(stderr, "Invalid bounds string, %s\n", s2);
        free(s);
        return;
    }
    for(i=1; i<4; i++) {
        p = strtok(NULL, ",");
        v = strtol(p, &end, 10);
        if((errno == ERANGE && (v == LONG_MAX || v == LONG_MIN)) || (errno != 0 && v == 0) || end == p) v = -1;
        if(v>=0) vals[4*mult+i] = v;
        else {
            fprintf(stderr, "Invalid bounds string, %s\n", s2);
            free(s);
            return;
        }
    }
    free(s);
}

void extract_usage() {
    fprintf(stderr, "\nUsage: PileOMeth extract [OPTIONS] reference.fa file.bam\n");
    fprintf(stderr,
"\n"
"Options:\n"
" -q INT           Minimum MAPQ threshold to include an alignment (default 10)\n"
" -p INT           Minimum Phred threshold to include a base (default 5). This\n"
"                  must be >0.\n"
" -D INT           Maximum per-base depth (default 2000)\n"
" -r STR           Region string in which to extract methylation\n"
" -l FILE          A BED file listing regions for inclusion. Note that unlike\n"
"                  samtools mpileup, this option will utilize the strand column\n"
"                  (column 6) if present. Thus, if a region has a '+' in this\n"
"                  column, then only metrics from the top strand will be\n"
"                  output. Note that the -r option can be used to limit the\n"
"                  regions of -l.\n"
" -o, --opref STR  Output filename prefix. CpG/CHG/CHH metrics will be\n"
"                  output to STR_CpG.bedGraph and so on.\n"
" --keepDupes      By default, any alignment marked as a duplicate is ignored.\n"
"                  This option causes them to be incorporated.\n"
" --keepSingleton  By default, if only one read in a pair aligns (a singleton)\n"
"                  then it's ignored.\n"
" --keepDiscordant By default, paired-end alignments with the properly-paired bit\n"
"                  unset in the FLAG field are ignored. Note that the definition\n"
"                  of concordant and discordant is based on your aligner\n"
"                  settings.\n"
" --noCpG          Do not output CpG methylation metrics\n"
" --CHG            Output CHG methylation metrics\n"
" --CHH            Output CHH methylation metrics\n"
" --OT INT,INT,INT,INT Inclusion bounds for methylation calls from reads/pairs\n"
"                  origination from the original top strand. Suggested values can\n"
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
"                  strands, respectively.\n");
}

int extract_main(int argc, char *argv[]) {
    char *opref = NULL, *oname, *p;
    int c, i;
    Config config;

    //Defaults
    config.keepCpG = 1; config.keepCHG = 0; config.keepCHH = 0;
    config.minMapq = 10; config.minPhred = 5; config.keepDupes = 0;
    config.keepSingleton = 0, config.keepDiscordant = 0;
    config.maxDepth = 2000;
    config.fai = NULL;
    config.fp = NULL;
    config.bai = NULL;
    config.reg = NULL;
    config.bedName = NULL;
    config.bed = NULL;
    for(i=0; i<16; i++) config.bounds[i] = 0;

    static struct option lopts[] = {
        {"opref",        1, NULL, 'o'},
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
        {"help",         0, NULL, 'h'}
    };
    while((c = getopt_long(argc, argv, "q:p:r:l:o:D:", lopts,NULL)) >= 0) {
        switch(c) {
        case 'h' :
            extract_usage();
            return 0;
        case 'o' :
            opref = strdup(optarg);
            break;
        case 'D' :
            config.maxDepth = atoi(optarg);
            break;
	case 'r':
	    config.reg = strdup(optarg);
	    break;
        case 'l' :
            config.bedName = optarg;
            break;
        case 1 :
            config.keepCpG = 0;
            break;
        case 2 :
            config.keepCHG = 1;
            break;
        case 3 :
            config.keepCHH = 1;
            break;
        case 4 :
            config.keepDupes = 1;
            break;
        case 5 :
            config.keepSingleton = 1;
            break;
        case 6 :
            config.keepDiscordant = 1;
            break;
        case 7 :
            parseBounds(optarg, config.bounds, 0);
            break;
        case 8 :
            parseBounds(optarg, config.bounds, 1);
            break;
        case 9 :
            parseBounds(optarg, config.bounds, 2);
            break;
        case 10 :
            parseBounds(optarg, config.bounds, 3);
            break;
        case 'q' :
            config.minMapq = atoi(optarg);
            break;
        case 'p' :
            config.minPhred = atoi(optarg);
            break;
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

    //Is there still a metric to output?
    if(!(config.keepCpG + config.keepCHG + config.keepCHH)) {
        fprintf(stderr, "You haven't specified any metrics to output!\nEither don't use the --noCpG option or specify --CHG and/or --CHH.\n");
        return -1;
    }

    //Open the files
    if((config.fai = fai_load(argv[optind])) == NULL) {
        fprintf(stderr, "Couldn't open the index for %s!\n", argv[optind]);
        extract_usage();
        return -2;
    }
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

    //Output files
    config.output_fp = malloc(sizeof(FILE *) * 3);
    if(opref == NULL) {
        opref = strdup(argv[optind+1]);
        p = strrchr(opref, '.');
        if(p != NULL) *p = '\0';
        fprintf(stderr, "writing to prefix:'%s'\n", opref);
    }
    oname = malloc(sizeof(char) * (strlen(opref)+14));
    if(config.keepCpG) {
        sprintf(oname, "%s_CpG.bedGraph", opref);
        config.output_fp[0] = fopen(oname, "w");
        if(config.output_fp[0] == NULL) {
            fprintf(stderr, "Couldn't open the output CpG metrics file for writing! Insufficient permissions?\n");
            return -3;
        }
        fprintf(config.output_fp[0], "track type=\"bedGraph\" description=\"%s CpG methylation levels\"\n", opref);
    }
    if(config.keepCHG) {
        sprintf(oname, "%s_CHG.bedGraph", opref);
        config.output_fp[1] = fopen(oname, "w");
        if(config.output_fp[1] == NULL) {
            fprintf(stderr, "Couldn't open the output CHG metrics file for writing! Insufficient permissions?\n");
            return -3;
        }
        fprintf(config.output_fp[1], "track type=\"bedGraph\" description=\"%s CHG methylation levels\"\n", opref);
    }
    if(config.keepCHH) {
        sprintf(oname, "%s_CHH.bedGraph", opref);
        config.output_fp[2] = fopen(oname, "w");
        if(config.output_fp[2] == NULL) {
            fprintf(stderr, "Couldn't open the output CHH metrics file for writing! Insufficient permissions?\n");
            return -3;
        }
        fprintf(config.output_fp[2], "track type=\"bedGraph\" description=\"%s CHH methylation levels\"\n", opref);
    }

    //Run the pileup
    extractCalls(&config);

    //Close things up
    hts_close(config.fp);
    fai_destroy(config.fai);
    if(config.keepCpG) fclose(config.output_fp[0]);
    if(config.keepCHG) fclose(config.output_fp[1]);
    if(config.keepCHH) fclose(config.output_fp[2]);
    hts_idx_destroy(config.bai);
    free(opref);
    if(config.reg) free(config.reg);
    if(config.bed) destroyBED(config.bed);
    free(oname);
    free(config.output_fp);

    return 0;
}
