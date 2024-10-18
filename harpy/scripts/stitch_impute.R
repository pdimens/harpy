suppressWarnings(suppressPackageStartupMessages(library("STITCH")))

# Params pulled in from Snakemake
bamlist <- snakemake@input[["bamlist"]]
posfile <- snakemake@input[["infile"]]
chr <- gsub(".stitch", "", basename(posfile))
outvcf <- snakemake@output[["vcf"]]
outdir <- normalizePath(dirname(outvcf))
outfile <- basename(outvcf)
logfile <- file(snakemake@log[["logfile"]], open = "wt")
# create tmpdir if not already present
tmpdr <- paste(outdir, "tmp", sep = "/")
dir.create(file.path(tmpdr), showWarnings = FALSE)

# model parameters
modeltype <- snakemake@params[["model"]]
K <- snakemake@params[["k"]]
S <- snakemake@params[["s"]]
bx <- as.logical(snakemake@params[["usebx"]])
bxlim <- snakemake@params[["bxlimit"]]
nGenerations <- snakemake@params[["ngen"]]
nCores <- snakemake@threads
inputBundleBlockSize <- NA
cli_args <- list(
    method               = modeltype,
    posfile              = posfile,
    bamlist              = bamlist,
    nCores               = nCores,
    nGen                 = nGenerations,
    chr                  = chr,
    K                    = K,
    S                    = S,
    use_bx_tag           = bx,
    bxTagUpperLimit      = bxlim,
    niterations          = 40,
    switchModelIteration = 39,
    splitReadIterations  = NA,
    outputdir            = outdir,
    output_filename      = outfile,
    tempdir              = tmpdr
)
# if there are any extra arguments provided to harpy by the -x argument
extra <- snakemake@params[["extra"]]
if(extra != ""){
    # convert the extra arguments into proper R types
    # converts numbers to numeric, vectors to vectors, leaves strings as-is
    interpret <- function(x){
    tryCatch({eval(parse(text=x))}, error = function(y){x})
    }
    extraargvals <- unlist(strsplit(extra, ","))
    n <- length(extraargvals)
    res <- c()
    startfrom <- 1
    for(i in 1:n){
        if(startfrom > i){
            next
        }
    currentstring <- extraargvals[i]
    if(grepl("c\\(", currentstring)){
        for(j in i:n){
            if(grepl("\\)", extraargvals[j])){
                .res <- paste(extraargvals[i:j], collapse = ",")
                startfrom <- j+1
                break
            }
        }
    } else {
        .res <- currentstring
    }
    res <- c(res, trimws(.res, "both"))
    }
    extra <- lapply(strsplit(res, "="), function(x){trimws(x, "both")})
    argnames <- sapply(extra,"[[",1)
    extraparams <- sapply(extra,"[[",2)
    extra_args <- lapply(extraparams, interpret)
    names(extra_args) <- argnames
    # combine the parameter list with the extra list, overwrite conflicts
    stitch_args <- utils::modifyList(cli_args, extra_args)
} else {
    stitch_args <- cli_args
}

sink(logfile, type = "output")
sink(logfile, type = "message")
do.call(STITCH, stitch_args)
sink()

debugdir <- paste(outdir, "debug", sep = "/")
# remove debug if it's empty
if (length(list.files(debugdir)) < 1) {
    unlink(debugdir, recursive = TRUE)
}