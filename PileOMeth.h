#include <inttypes.h>
#include <stdio.h>
#include <zlib.h>
#include "htslib/sam.h"
#include "htslib/faidx.h"

/*! @typedef
 @abstract Structure to hold one region defined in a BED file
 @field	tid	The chromosome ID, defined by bam_hdr_t
 @field start	0-based start position
 @field end	1-based end position
 @field	strand	0: Ignore strand (".")
                1: Top strand ("+")
                2: Bottom strand ("-")

 @discussion The start and end coordinates will be truncated if they go beyond
 the bounds of a chromosome, as indicated by bam_hdr_t.
*/
typedef struct {
    int32_t tid, start, end; //0-based, [start, end)
    int8_t strand; //0: ., 1: +, 2: -
} bedRegion;

/*! @typedef
 @abstract Structure to hold one region defined in a BED file
 @field region	Pointer to the regions
 @field	n	Number of regions
 @field	m	maxmum number of regions that can currently be stored.

 @discussion You must free this with destroyBED()
*/
typedef struct {
    bedRegion *region;
    int32_t n, m; //Current (n) and maximal possible (m) number of regions
} bedRegions;

/*! @typedef
 @abstract Global configuration information structure
 @field	keepCpG	1: Output CpG metrics
 @field keepCHG	1: Output CHG metrics
 @field	keepCHH	1: Output CHH metrics
 @field	minMapq	Minimum MAPQ value to include a read (-q)
 @field	minPhred	Minimum Phred score to include a base (-p)
 @field	keepDupes	1: Include marked duplicates when calculating metrics
 @field	maxDepth	Maximum depth for the pileup
 @field noDiscordant	1: Do not include discordantly aligned reads when calculating metrics
 @field	noSingleton	1: Do not include singletons when calculating metrics
 @field output_fp	Output file pointers (to CpG, CHG, and CHH, respectively)
 @field	reg	A region specified by -r
 @field fp	Input file pointer (must be a BAM or CRAM file)
 @field	bai	The index for fp
 @field bed	Pointer to regions specified in a BED file (-l option)
 @field fai	Fasta file index pointer
*/
typedef struct {
    int keepCpG, keepCHG, keepCHH;
    int minMapq, minPhred, keepDupes, maxDepth;
    int noDiscordant, noSingleton;
    FILE **output_fp;
    char *reg;
    htsFile *fp;
    hts_idx_t *bai;
    char *bedName;
    bedRegions *bed;
    faidx_t *fai;
} Config;

/*! @function
 @abstract Determine the strand from which a read arose
 @param	b	Pointer to an alignment
 @returns	1, original top; 2, original bottom; 3, complementary to the original top; 4, complementary to the original bottom
 @discussion There are two methods used to determine the strand of origin. The 
 first method uses XG auxiliary tags as output by Bison and Bismark. For details
 on how strand is determined from this, see the source code for this function or
 the documentation for either Bison or bismark. This method can handle non-
 directional libraries. If an XG auxiliary tag is not present, then the
 orientation of each alignment (and whethers it's read #1 or #2 of a pair, if
 applicable) determines the strand.
*/
int getStrand(bam1_t *b);

int posOverlapsBED(int32_t tid, int32_t pos, bedRegions *regions, int idxBED);
int spanOverlapsBED(int32_t tid, int32_t start, int32_t end, bedRegions *regions, int *idx);
int readStrandOverlapsBED(bam1_t *b, bedRegion region);
void sortBED(bedRegions *regions);
void destroyBED(bedRegions *regions);
bedRegions *parseBED(char *fn, bam_hdr_t *hdr);
