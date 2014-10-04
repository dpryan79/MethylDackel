PileOMeth (a temporary name derived due to it using a PILEup to extract METHylation metrics) will process a coordinate-sorted and indexed BAM or CRAM file containing some form of BS-seq alignments and extract per-base methylation metrics from them. PileOMeth requires an indexed fasta file containing the reference genome as well.

Prerequisites
=============

PileOMeth requires HTSlib version 1.1 or newer (version 1.0 is currently unsupported).

Compilation
===========

PileOMeth can either dynamically link to HTSlib or statically incorporate it. To compile for dynamic linkage, use:

`gcc -Wall -O3 -I/path/to/headers -L/path/to/libhts.so -o PileOMeth PileOMeth.c -lhts -lz -lpthread`

If your HTSlib headers are under ~/include/htslib then use `-I~/include`.

If you prefer to statically incorporate HTSlib (e.g., due to this being the only program needing it), compile HTSlib and then use the following:

`gcc -Wall -O3 -I/path/to/htslib/compilation -o PileOMeth PileOMeth.c /path/to/htslib/compilation/libhts.a -lz -lpthread`

If you downloaded and compiled HTSlib in ~/Downloads/htslib-1.1, then `/path/to/htslib/compilation` is `~/Downloads/htslib-1.1`.

Usage
=====

The most basic usage of PileOMeth is as follows:

`PileOMeth reference_genome.fa alignments.bam`

This will calculate per-base CpG metrics and output them to `alignments_CpG.bedGraph`, which is a standard bedGraph file with column 4 being the number of reads/read pairs with evidence for a methylated C at a given position and column 5 the equivalent for an unmethylated C. An alternate output filename prefix can be specified with the `-o some_new_prefix` option.

By default, PileOMeth will only calculate CpG metrics, but CHG and CHH metrics are supported as well (see the --CHH and --CHG options). If you would like to ignore CpG, metrics, simply specify --noCpG. Each type of metric is output to a different file.

PileOMeth can filter reads and bases according to MAPQ and Phred score, respectively. The default minimums are MAPQ >= 10 and Phred >= 5, though these can be changed with the -q and -p options.

To do list
==========

 * Add an -l option to parse only regions in a BED file
 * add an -r/--region REGION option to just look at one region
 * Enable trimming based on M-bias
 * The -D option seems to only be approximate. Is the an htslib issue?
 * Is it possible to support non-directional libraries with this method?
 * If a BAM (or ideally CRAM) file isn't yet indexed, we should do that automatically
 * Is the output format the most convenient (it's what Bison uses, so converters have already been written)? It makes more sense to output a predefined VCF format, which would allow processing multiple samples at once. This would require a spec., which should have pretty broad input.
 * Perhaps things should be restructured to allow making a library out of this, for easier incorporation into python (e.g., pysam).
 * Test to ensure that the results are correct!
