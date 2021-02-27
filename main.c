#include <string.h>
#include <stdio.h>
#include "htslib/hts.h"
#include "version.h" //This has the VERSION define

pthread_mutex_t positionMutex;
pthread_mutex_t bwMutex;
pthread_mutex_t outputMutex;
uint32_t globalTid = 0;
uint32_t globalPos = 0;
uint32_t globalEnd = 0;
uint32_t bin = 0;
uint32_t outputBin = 0;
uint64_t globalnVariantPositions = 0;

int mbias_main(int argc, char *argv[]);
int extract_main(int argc, char *argv[]);
int mergeContext_main(int argc, char *argv[]);
int perRead_main(int argc, char *argv[]);
void print_version(void);

void usage_main() {
    fprintf(stderr, "MethylDackel: A tool for processing bisulfite sequencing alignments.\n"
"Version: %s (using HTSlib version %s)\n", VERSION, hts_version());
    fprintf(stderr,
"Usage: MethylDackel <command> [options]\n\n"
"Commands:\n"
"    mbias    Determine the position-dependent methylation bias in a dataset,\n"
"             producing diagnostic SVG images.\n"
"    extract  Extract methylation metrics from an alignment file in BAM/CRAM\n"
"             format.\n"
"    mergeContext   Combine single Cytosine metrics from 'MethylDackel extract' into\n"
"             per-CpG/CHG metrics.\n"
"    perRead  Generate a per-read methylation summary.\n"
);
}

int main(int argc, char *argv[]) {
    if(argc == 1) {
        usage_main();
        return 0;
    } else if(strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        usage_main();
        return 0;
    } else if(strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "--version") == 0) {
        print_version();
        return 0;
    } else if(strcmp(argv[1], "extract") == 0) {
        return extract_main(argc-1, argv+1);
    } else if(strcmp(argv[1], "mbias") == 0) {
        return mbias_main(argc-1, argv+1);
    } else if(strcmp(argv[1], "mergeContext") == 0) {
        return mergeContext_main(argc-1, argv+1);
    } else if(strcmp(argv[1], "perRead") == 0) {
        return perRead_main(argc-1, argv+1);
    } else {
        fprintf(stderr, "Unknown command!\n");
        usage_main();
        return -1;
    }
}
