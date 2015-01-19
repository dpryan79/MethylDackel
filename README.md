PileOMeth (a temporary name derived due to it using a PILEup to extract METHylation metrics) will process a coordinate-sorted and indexed BAM or CRAM file containing some form of BS-seq alignments and extract per-base methylation metrics from them. PileOMeth requires an indexed fasta file containing the reference genome as well.

Prerequisites
=============

PileOMeth requires HTSlib version 1.1 or newer (version 1.0 is currently unsupported). As HTSlib is a submodule of this repository, this prerequisite should be automatically dealt with.

Compilation
===========

Compilation and installation can be performed via:

    make
    make install path=/some/installation/path

As HTSlib is now a submodule of this repository, you no longer need to manually download and compile it.

Usage
=====

The most basic usage of PileOMeth is as follows:

`PileOMeth reference_genome.fa alignments.bam`

This will calculate per-base CpG metrics and output them to `alignments_CpG.bedGraph`, which is a standard bedGraph file with column 4 being the number of reads/read pairs with evidence for a methylated C at a given position and column 5 the equivalent for an unmethylated C. An alternate output filename prefix can be specified with the `-o some_new_prefix` option.

By default, PileOMeth will only calculate CpG metrics, but CHG and CHH metrics are supported as well (see the --CHH and --CHG options). If you would like to ignore CpG, metrics, simply specify --noCpG. Each type of metric is output to a different file.

PileOMeth can filter reads and bases according to MAPQ and Phred score, respectively. The default minimums are MAPQ >= 10 and Phred >= 5, though these can be changed with the -q and -p options.

To do list
==========

 * Enable trimming based on M-bias
 * Is the output format the most convenient (it's what Bison uses, so converters have already been written)? It makes more sense to output a predefined VCF format, which would allow processing multiple samples at once. This would require a spec., which should have pretty broad input.
 * Need to finish restructuring things to allow easy library incorporation
