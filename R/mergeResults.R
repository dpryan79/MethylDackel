library(rtracklayer)

## ideally, we'd like to do something along the lines of:
##
##   res <- mergeResults(BAMfiles, mode="union", minCounts=5, how="perGroup") 
##
## or something along those lines.  Then the same code in Python or R would
## do the same thing by wrapping PileOMeth and running through the same list.
##
BAMfiles <- list.files(patt=".bam$")

## later on, we'll want to generate these in-flight from BAMs, using PileOMeth
## for now, just open them up 
##
options("mc.cores"=4)
getMethAndCounts <- function(bamFiles, minSamples=NULL, exptName="tabulated"){
  library(parallel)
  if (is.null(minSamples)) minSamples <- length(bamFiles)
  countFiles <- sub(".bam$", "_CpG.counts.bw", bamFiles)
  names(countFiles) <- sapply(countFiles, function(x) strsplit(x,"\\.")[[1]][1])
  names(countFiles) <- gsub("-", "_", names(countFiles))
  countRanges <- GRangesList(mclapply(countFiles, import))
  coords <- sort(unname(unique(unlist(countRanges))))
  coords <- coords[ countOverlaps(coords, countRanges) >= minSamples ]
  counts <- do.call(cbind,
                    lapply(countRanges, 
                           function(cr) subsetByOverlaps(cr, coords)$score))
  rownames(counts) <- paste0(seqnames(coords), ":", start(coords))
  methFiles <- sub("_CpG.counts.bw", "_CpG.meth.bw", countFiles)
  methRanges <- GRangesList(mclapply(methFiles, import, selection=coords))
  meths <- do.call(cbind, lapply(methRanges, function(mr) mr$score))
  rownames(meths) <- paste0(seqnames(coords), ":", start(coords))

  ## reorder to match input file order...
  counts <- counts[ , names(countFiles) ]
  meths <- meths[ , names(countFiles) ]

  stopifnot(identical(colnames(meths), colnames(counts)))
  write.table(meths, paste(exptName, "methylation", "txt", sep="."))
  write.table(counts, paste(exptName, "counts", "txt", sep="."))
  invisible(list(counts=counts, methylation=meths))
}
