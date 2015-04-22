#include <string.h>
#include <stdio.h>
#include "htslib/hts.h"
#include "version.h" //This has the VERSION define

int mbias_main(int argc, char *argv[]);
int extract_main(int argc, char *argv[]);
int extract_fraction(int argc, char *argv[]);
int extract_counts(int argc, char *argv[]);
int extract_logit(int argc, char *argv[]);

void usage_main() {
    fprintf(stderr, "PileOMeth: A tool for processing bisulfite sequencing alignments.\n"
"Version: %s (using HTSlib version %s)\n", VERSION, hts_version());
    fprintf(stderr,
"Usage: PileOMeth <command> [options]\n\n"
"Commands:\n"
"    mbias    Determine the position-dependent methylation bias in a dataset,\n"
"             producing diagnostic SVG images.\n"
"    extract  Extract methylation metrics from an alignment file in BAM/CRAM\n"
"             format.\n"
);
}

int main(int argc, char *argv[]) {
    if(argc == 1) {
        usage_main();
        return 0;
    } else if(strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        usage_main();
        return 0;
    } else if(strcmp(argv[1], "extract") == 0) {
        return extract_main(argc-1, argv+1);
    } else if(strcmp(argv[1], "mbias") == 0) {
        return mbias_main(argc-1, argv+1);
    } else {
        fprintf(stderr, "Unknown command!\n");
        usage_main();
        return -1;
    }
}
