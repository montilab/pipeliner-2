#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Biobase))

args <- commandArgs(trailingOnly=TRUE)
data <- read.delim(args[[1]], sep="\t")

fdat <- data.frame(id=data[,1])
rownames(fdat) <- fdat$id

pdat <- data.frame(id=colnames(data)[2:ncol(data)])
rownames(pdat) <- pdat$id

edat <- data[,2:ncol(data)]
rownames(edat) <- rownames(fdat)
edat <- as.matrix(edat)

fdata <- new("AnnotatedDataFrame", data=fdat)
pdata <- new("AnnotatedDataFrame", data=pdat)
eset <- ExpressionSet(assayData=edat, phenoData=pdata, featureData=fdata)

saveRDS(eset, "eset.rds")