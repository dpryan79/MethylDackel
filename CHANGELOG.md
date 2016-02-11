Version 0.1.10:

   * Added the --version/-v option.
   * Fixed the -h/--help option, which actually wasn't recognized before!
   * Added the --minDepth option.

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
