#!/usr/bin/env Rscript

featurecounts_importer = function(fnames) {
  first <- read.delim(fnames[1], header=TRUE, sep="\t", skip=1, stringsAsFactors=FALSE)
  genes <- first[[1]]
  rest <- do.call(cbind, lapply(fnames, function(fname){
    df <- read.delim(fname, header=TRUE, sep="\t", skip=1, stringsAsFactors=FALSE)
    ## trying to make it more robust
    count_idx <- grep(".bam", colnames(df))
    stopifnot(length(count_idx)==1)
    stopifnot(count_idx>0)
    df[,count_idx]
    ## old version
    #df[,ncol(df)]
  }))
  colnames(rest) <- basename(fnames)
  rest <- data.frame(genes, rest, stringsAsFactors=FALSE)
  return(rest)
}

# Parse arguments
args <- commandArgs(trailingOnly=TRUE)
stopifnot(args[[1]] == "--format")
format <- as.character(args[[2]])
files <- args[3:length(args)]

if (format == "featurecounts") {
    mat <- featurecounts_importer(files)
    write.table(mat, file="counts.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}
