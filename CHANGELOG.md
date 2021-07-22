Version 0.6.0:

   * Added the `minConversionEfficiency` option, which allows filtering out incompletely converted alignments on the fly. Note that doing this is generally NOT recommended. (issue #61)

Version 0.5.3:

   * Fixed an issue with the `perRead` subcommand, wherein the requireFlags option didn't fully work (a read would pass if it had at least one of the required flags set, rather than all of them). (issue #117)

Version 0.5.2:

   * Rewrote how read-pair overlap handling is performed. It now uses the constructor/destructor mechanism from htslib instead of using internal htslib structures and functions. This allows supporting newer htslib versions. Currently 1.11 is the only tested and working version, due to changes in the pileup constructor interface in it. (issue #99)

Version 0.5.1:

   * Fixed an issue in `MethylDackel mbias` due to an uninitialized value (issue #93).

Version 0.5.0:

   * Fixed an issue with the `--cytosine_report` option where the reported chromosome name could be wrong IF the input BAM files were very sparse and multiple threads were used. (issue #88)
   * libBigWig is now an external dependency as it's needed for handling the "binary bismap format"
   * Added support for blacklisting reads according to a "binary bismap file". See [here](https://github.com/dpryan79/MethylDackel/blob/master/BBM_Specification.md) for details. Code for this was contributed by @valiec and @bwlang at New England Biolabs.

Version 0.4.0:

   * Switched to an external htslib. It is currently compatible with htslib versions 1.4 through 1.9 (the latest one at the time of release).
   * Added the `--perRead` option, as used in https://www.biorxiv.org/content/10.1101/481283v1
   * Removed the `--maxDepth` option, it's not longer required.
   * Fixed issue #58, the `--keepDupes` flag now changes `--ignoreFlags`.
   * Fixed issue #59, the confidence intervals should no longer extend outside of [0, 1].

Version 0.3.0:

   * Added `--nOT`, `--nOB`, `--nCTOT` and `--nCTOB`, which are related to `--OT` and similar. The difference is that the specified number of bases will be ignored, regardless of read length. So `--nOT 5,5,5,5` will result in the 5 bases on either end of both reads being ignored.

Version 0.2.1:

   * The `--methylKit` option output methylation values between 0-1 rather than 0-100. This has been corrected.

Version 0.2.0:

   * Changed the package name from PileOMeth to MethylDackel. It is unfortunate that the temporary "PileOMeth" name stuck around for so long.
   * Fixed the plotting where the sometimes the plotted lines escape the bounds of the graph. The cause for this was that read #2 was being ignored when the graph bounds were being computed.
   * Added the `--requireFlags`/`-r` option, which is equivalent to the -f option in samtools. The default is 0, which requires nothing.

Version 0.1.13:

   * Added the `--ignoreFlags`/`-F` option, which is equivalent to the -F option in samtools. The default is 0xF00, which ignores duplicates, QC fail, secondary, and supplemental alignments.

Version 0.1.12:

   * Fixed handling of hard-clipped bases, which caused a segfault before.

Version 0.1.11:

   * The --minDepth/-d option is actually recognized now. Sorry about that!

Version 0.1.10:

   * Added the --version/-v option.
   * Fixed the -h/--help option, which actually wasn't recognized before!
   * Added the --minDepth option.
   * Fixed handling of the `XG` auxiliary tag, since some aligners other than bismark are using it (and for different purposes). This is the single biggest annoyance of custom auxiliary tags.

Version 0.1.9:

   * Added the --methylKit option
   * Added automated testing with `make test`

Version 0.1.8:

   * Fix compilation issues with clang

Version 0.1.7:

   * BED files are now sorted properly before extraction, so nothing is skipped now.
   * The entire BED file isn't iterated over when there's no overlap. This significantly speed things up.
   * The first alignment overlapping a given BED interval is no longer ignored.

Version 0.1.6 :

   * Fixed a bug in the --mergeContext option that sometimes resulted in incorrect mergin and duplicate lines.

Version 0.1.5 :

    * htslib is no longer a submodule, but instead version 1.2.1 is directly included. This enables easier integration into Galaxy.

Version 0.1.4 :

    * Fixed a bug in getStrand() that handled SE alignments incorrectly.

Version 0.1.3 :

    * Added the --mergeContext option to `PileOMeth extract`.
    * Added the `PileOMeth mergeContext` command. This can only be used for the default output.
    * Added this file.
    * The header line of each output file is now modified according to the presence of --logit, --fraction, and --counts.

Version 0.1.2 :

    * Fixed a mistake in the Makefile, where "prefix" was used rather than "path".
    * Added the --logit, --counts, and --fraction option, thanks to Tim Triche

Version 0.1.1 :

    * Document the output format in README.md.
    * Convert the 4th column in the output file(s) to contain the methylation level from 0-100, rather than 0-1000.

Version 0.1.0 :

    * First actual release. The code seems to be stable enough now to be usable.
