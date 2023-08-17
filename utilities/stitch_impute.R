#!/usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(library("STITCH")))

# Params pulled in from Snakemake
bamlist <- snakemake@input[["bamlist"]]
posfile <- snakemake@input[["infile"]]
chr <- gsub(".stitch", "", basename(posfile))
outdir <- normalizePath(dirname(snakemake@output[[1]]))
outfile <- basename(snakemake@output[[1]])
logfile <- file(snakemake@log[[1]], open = "wt")
# create tmpdir if not already present
tmpdr <- paste(outdir, "tmp", sep = "/")
dir.create(file.path(tmpdr), showWarnings = FALSE)

# model parameters
parameters <- snakemake@params[["parameters"]]
#extra <- snakemake@params[["extra"]]
modeltype <- parameters$model
K <- parameters$k
S <- parameters$s
.bx <- toupper(parameters$usebx)
bx <- .bx == "TRUE" || .bx == "YES" || .bx == "Y"
bxlim <- parameters$bxlimit
nGenerations <- parameters$ngen
nCores <- snakemake@threads
inputBundleBlockSize <- NA

# WTF is a genfile?
sink(logfile, type = "output")
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
    bxTagUpperLimit         = bxlim,
    niterations             = 40,
    switchModelIteration    = 39,
    splitReadIterations     = NA,
    outputdir               = outdir,
    output_filename         = outfile,
    tempdir                 = tmpdr
)
sink()

debugdir <- paste(outdir, "debug", sep = "/")
# remove debug if it's empty
if (length(list.files(debugdir)) < 1) {
    unlink(debugdir, recursive = TRUE)
}