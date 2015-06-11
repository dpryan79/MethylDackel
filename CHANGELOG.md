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
