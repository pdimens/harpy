#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("STITCH"))


# Params pulled in from Snakemake
bamlist <- snakemake@input[["bamlist"]]
chr <- readLines(snakemake@input[["chromosome"]])[1]
posfile <- snakemake@input[["infile"]]
outdir <- normalizePath(dirname(snakemake@output[[1]]))
outfile <- basename(snakemake@output[[1]])
logfile <- file(snakemake@log[[1]], open = "wt")

# model parameters 
parameters <- snakemake@params[["parameters"]]
modeltype <- parameters$model
K <- parameters$k
S <- parameters$s
.bx <- toupper(parameters$useBX) 
bx <- .bx == "TRUE" || .bx == "YES" || .bx == "Y"
nGenerations <- parameters$nGen
nCores <- snakemake@threads
inputBundleBlockSize <- NA

# WTF is a genfile?
sink(logfile ,type = "output")
sink(logfile, type = "message")
STITCH(
    method                  = modeltype,
    posfile                 = posfile,
    bamlist                 = bamlist,
    nCores                  = nCores,
    nGen                    = nGenerations,
    chr                     = chr,
    K                       = K,
    S                       = S,
    use_bx_tag              = bx,
    bxTagUpperLimit         = 50000,
    niterations             = 40,
    switchModelIteration    = 39,
    splitReadIterations     = NA,
    outputdir               = outdir,
    output_filename         = outfile
)
sink()
sink()
