[![Master build status](https://travis-ci.org/dpryan79/PileOMeth.svg?branch=master)](https://travis-ci.org/dpryan79/PileOMeth)
[![Galaxy installation](https://img.shields.io/badge/install%20with-galaxy-brightgreen.svg?style=flat-square)](https://toolshed.g2.bx.psu.edu/)

PileOMeth (a temporary name derived due to it using a PILEup to extract METHylation metrics) will process a coordinate-sorted and indexed BAM or CRAM file containing some form of BS-seq alignments and extract per-base methylation metrics from them. PileOMeth requires an indexed fasta file containing the reference genome as well.

Prerequisites
=============

PileOMeth now includes HTSlib version 1.2.1, though any version starting with 1.1 should work. This was originally a submodule, but that made things difficult for Galaxy integration.

Compilation
===========

Compilation and installation can be performed via:

    git clone https://github.com/dpryan79/PileOMeth.git
    cd PileOMeth
    make
    make install prefix=/some/installation/path


Methylation Context
===================

PileOMeth groups all Cytosines into one of three sequence contexts: [CpG](http://en.wikipedia.org/wiki/CpG_site), CHG, and CHH. Here, H is the IUPAC ambiguity code for any nucleotide other than G. If an N is encountered in the reference sequence, then the context will be assigned to CHG or CHH, as appropriate (e.g., CNG would be categorized as in a CHG context and CNC as in a CHH context). If a Cytosine is close enough to the end of a chromosome/contig such that its context can't be inferred, then it is categorized as CHH (e.g., a Cytosine as the last base of a chromosome is considered as being in a CHH context).

Usage
=====

The most basic usage of PileOMeth is as follows:

    PileOMeth extract reference_genome.fa alignments.bam

This will calculate per-base CpG metrics and output them to `alignments_CpG.bedGraph`, which is a standard bedGraph file with column 4 being the number of reads/read pairs with evidence for a methylated C at a given position and column 5 the equivalent for an unmethylated C. An alternate output filename prefix can be specified with the `-o some_new_prefix` option.

By default, PileOMeth will only calculate metrics for Cytosines in a CpG context, but metrics for those in CHG and CHH contexts are supported as well (see the --CHH and --CHG options). If you would like to ignore Cytosines in CpGs, simply specify --noCpG. Each type of metric is output to a different file. For per-CpG and/or per-CHG rather than per-Cytosine metrics, see "Per-CpG/CHG metrics", below.

PileOMeth can filter reads and bases according to MAPQ and Phred score, respectively. The default minimums are MAPQ >= 10 and Phred >= 5, though these can be changed with the -q and -p options. PileOMeth can also account for methylation bias (described below) with the `--OT`, `--OB`, `--CTOT`, and `--CTOB` options.

Single Cytosine methylation metrics extraction
==============================================

`PileOMeth extract` produces a variant of [bedGraph](http://genome.ucsc.edu/goldenpath/help/bedgraph.html) that's similar to the "coverage" file produced by [Bismark](http://www.bioinformatics.babraham.ac.uk/projects/bismark/) and [Bison](https://github.com/dpryan79/bison). In short, each line consists of 6 tab separated columns:

1. The chromosome/contig/scaffold name
2. The start coordinate
3. The end coordinate
4. The methylation percentage rounded to an integer
5. The number of alignments/pairs reporting methylated bases
6. The number of alignments/pairs reporting unmethylated bases

All coordinates are 0-based half open, which conforms to the bedGraph definition. When paired-end reads are aligned, it can often occur that their alignments overlap. In such cases, PileOMeth will not count both reads of the pair in its output, as doing so would lead to incorrect downstream statistical results.

An example of the output is below:

    track type="bedGraph" description="SRR1182519.sorted CpG methylation levels"
    1	25115	25116	100	3	0
    1	29336	29337	50	1	1

Note the header line, which starts with "track". The "description" field is used as a label in programs such as [IGV](http://www.broadinstitute.org/igv/). Each of the subsequent lines describe single Cytosines, the 25116th and 29337th base on chromosome 1, respectively. The first position has 3 alignments (or pairs of alignments) indicating methylation and 0 indicating unmethylation (100% methylation) and the second position has 1 alignment each supporting methylation and unmethylation (50% methylation).

Per-CpG/CHG metrics
===================

In many circumstances, it's desireable for metrics from individual Cytosines in a CpG to be merged, producing per-CpG metrics rather than per-Cytosine metrics. This can be accomplished with the `--mergeContext` option in `PileOMeth extract`. If this is used, then this output:

    track type="bedGraph" description="SRR1182519.sorted CpG methylation levels"
    1	25114	25115	100	2	1
    1	25115	25116	100	3	0

is changed to this:

    track type="bedGraph" description="SRR1182519.sorted merged CpG methylation levels"
    1	25114	25116	100	5	1

This also works for CHG-level metrics. If bedGraph files containing per-Cytosine metrics already exist, they can be converted to instead contain per-CpG/CHG metrics with `PileOMeth mergeContext`.

Excluding low-coverage regions
==============================

If your downstream analysis requires an absolute minimum coverage (here, defined as the number of methylation calls kept after filtering for MAPQ, phred score, etc.), you can use the `--minDepth` option to achieve this. By default, `PileOMeth extract` will output all methylation metrics as long as the coverage is at least 1. If you use `--minDepth 10`, then only sites covered at least 10x will be output. This works in conjunction with the `--mergeContext` option, above. So if you request per-CpG context output (i.e., with `--mergeContext`) and `--minDepth 10` then only CpGs with a minimum coverage of 10 will be output.

Logit, fraction, and counts only output
=======================================

The standard output described above can be modified if you supply the `--fraction`, `--counts`, or `--logit` options to `PileOMeth extract`.

The `--fraction` option essentially produces the first 4 columns of the standard output described above. The only other difference is that the range of the 4th column is now between 0 and 1, instead of 0 and 100. Instead of producing a file ending simply in `.bedGraph`, one ending in `.meth.bedGraph` will instead be produced.

The `--counts` option produces the first three columns of the standard output followed by a column of total coverage counts. This last column is equivalent to the sum of the 5th and 6th columns of the standard output. The resulting file ends in `.counts.bedGraph` rather than simply `.bedGraph`.

The `--logit` option produces the first three columns of the standard output followed by the logit transformed methylation fraction. The logit transformation is log(Methylation fraction/(1-Methylation fraction)). Note that log uses base e. Logit transformed methylation values range between +/- infinity, rather than [0,1]. The resulting file ends in `.logit.bedGraph` rather than simply `.bedGraph`.

Note that these options may be combined with `--mergeContext`. However, `PileOMeth mergeContext` can not be used after the fact to combine these.

methylKit-compatible output
===========================

methylKit has its own format, which can be produced with the `--methylKit` option. Merging Cs into CpGs or CHGs is forbidden in this format. Likewise, this option is mutually exclusive with `--logit` et al.

Methylation bias plotting and correction
========================================

In an ideal experiment, we expect that the probability of observing a methylated C is constant across the length of any given read. In practice, however, there are often increases/decreases in observed methylation rate at the ends of reads and/or more global changes. These are termed methylation bias and including such regions in the extracted methylation metrics will result in noisier and less accurate data. For this reason, users are strongly encouraged to make a methylation bias plot. PileOMeth comes with a function for just this purpose:

    PileOMeth mbias reference_genome.fa alignments.sorted.bam output_prefix

That command will create a methylation bias (mbias for short) plot for each of the strands for which there are valid alignments. The command can take almost all of the same options as `PileOMeth extract`, so if you're interested in looking at only a single region or only CHH and CHG metrics then you can do that (run `PileOMeth mbias -h` for the full list of options). The resulting mbias graphs are in SVG format and can be viewed in most modern web browsers:

![An example SVG image](https://rawgit.com/dpryan79/PileOMeth/master/example_OT.svg)

If you have paired-end data, both reads in the pair will be shown separately, as is the case above. The program will suggest regions for inclusion ("--OT 2,0,0,98" above) and mark them on the plot, if applicable. The format of this output is described in `PileOMeth extract -h`. These suggestions should not be accepted blindly; users are strongly encouraged to have a look for themselves and tweak the actual bounds as appropriate. The lines indicate the average methylation percentage at a given position and the shaded regions the 99.9% confidence interval around it. This is useful in gauging how many methylation calls a given position has relative to its neighbors. Note the spike in methylation at the end of read #2 and the corresponding dip at the beginning of read #1. This is common and these regions can be ignored with the suggested trimming bounds. Note also that the numbers refer to the first and last base that should be included during methylation extraction, not the last and first base to ignore!.

To do list
==========

 - [ ] Is the output format the most convenient (it's what Bison uses, so converters have already been written)? It makes more sense to output a predefined VCF format, which would allow processing multiple samples at once. This would require a spec., which should have pretty broad input.
