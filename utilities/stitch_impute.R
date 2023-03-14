#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("STITCH"))

# Params pulled in from Snakemake
bamlist <- snakemake@input[["bamlist"]]
posfile <- snakemake@input[["infile"]]
chr <- gsub(".stitch", "", basename(posfile))
outdir <- normalizePath(dirname(snakemake@output[[1]]))
outfile <- basename(snakemake@output[[1]])
logfile <- file(snakemake@log[[1]], open = "wt")

# model parameters 
parameters <- snakemake@params[["parameters"]]
#extra <- snakemake@params[["extra"]]
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

# Remove all the extra RData files and directories
unlink(paste(outdir, "input", sep = "/"), recursive = TRUE)
unlink(paste(outdir, "RData", sep = "/"), recursive = TRUE)
debugdir <- paste(outdir, "debug", sep = "/")
# remove debug if it's empty
if (length(list.files(debugdir)) < 1){
    unlink(debugdir, recursive = TRUE)
}