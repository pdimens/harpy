#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("STITCH"))


# Params pulled in from Snakemake
bamlist <- snakemake@input[["bamlist"]]
chr <- snakemake@input[["chromosome"]]
outdir <- normalizePath(dirname(snakemake@output[[1]]))
outfile <- basename(snakemake@output[[1]])
tempdir <- paste(outdir, ".tempdir", sep = "/")
# model parameters 
modeltype <- snakemake@params[["model"]]
K <- snakemake@params[["K"]]
S <- snakemake@params[["S"]]
bx <- snakemake@params[["useBarcodes"]] == "TRUE"
nGenerations <- snakemake@params[["nGenerations"]]
nCores <- snakemake@threads
inputBundleBlockSize <- NA

# WTF is a genfile?
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
    tempdir                 = tempdir,
    outputdir               = outdir,
    output_filename         = outfile
)

#q()
#
## Libraries
#library("STITCH")
#library("parallel")
#library("Rcpp")
#library("readr")
#library("tidyverse")
#
## Set directories for testing
#datadir <- file.path(dirpath, "15-imputation/")
#resultsdir <- file.path(dirpath, "15-imputation/test-results/")
#tempdir <- file.path(dirpath, "15-imputation/test-temp/")
#
#
## bamlist
#shad_bamlist <- paste0(datadir, "bamlist_chr1-3.txt")
#
## chromosome list
#shad_chr <- readr::read_delim(paste0(datadir, "chromosome_list.bed"), col_names = FALSE, delim=" ")
#chr.list<-as.list(shad_chr[[1]])
#
## STITCH parallel
#for (i in 1 : length(chr.list)) {
#  chr <-  chr.list[i]
#  outputdir[i]=paste(c(resultsdir,'ShadHap_',chr,'_K',shad_K,'_nGen',shad_nGen,'_S',shad_S,'/'),collapse='')
#  shad_posfile[i] <- paste0(datadir, 'position_files/', 'pos_file_',chr,'.txt')
# 
#  STITCH(method = "pseudoHaploid",
#         posfile = shad_posfile[i],
#         bamlist = shad_bamlist,
#         nCores = n_cores,
#         nGen = shad_nGen,
#         chr = chr,
#         K = shad_K,
#         S = shad_S,
#         use_bx_tag = TRUE,
#         bxTagUpperLimit = 50000,
#         niterations = 40,
#         switchModelIteration = 39,
#         splitReadIterations = NA,
#         tempdir = tempdir,
#         outputdir = outputdir[i])
#  }