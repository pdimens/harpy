import os
import re
import sys
from snakemake.utils import Paramspace
import pandas as pd

infile = config["infile"]
samplefile = config["samplefile"]
paramspace  = Paramspace(pd.read_csv(samplefile, sep="\t"), param_sep = "", filename_params="*")
#outprefix = "DATA"
bn = os.path.basename(infile)
inprefix = re.sub(r"_[IR][12].fastq.gz", "", bn)
inprefixfull = re.sub(r"_[IR][12].fastq.gz", "", infile)
infiles = [f"{prefixfull}_{i}.fastq.gz" for i in ["I1", "I2","R1","R2"]]
indir = os.path.dirname(infile)
outdir = f"Demultiplex/{inprefix}/"

# mv functions to harpy executable?
def checkfiles(prefix, prefixfull):
    filelist = []
    report = False
    for i in ["I1", "I2","R1","R2"]:
        chkfile = f"{prefixfull}_{i}.fastq.gz"
        TF = os.path.exists(chkfile)
        report = True if not TF else report
        symbol = " " if TF else "X"
        filelist.append(f"\033[91m{symbol}\033[0m  {prefix}_{i}.fastq.gz")
    if report:
        print(f"Not all necessary files with prefix {prefix} present:")
        _ = [print(i, file = sys.stderr) for i in filelist]
        # EXIT WITH ERROR HERE


def get_samplenames(smpl):
    with open(smpl, "r") as f:
        rows = [i.split("\t")[0] for i in f.readlines()]
        # rm row of column names
        rows.pop(0)
        return rows


rule link_reads:
    input:
        indir + "/" + inprefix + "{part}" + ".fastq.gz"
    output:
        temp(outdir + "DATA_{part}.fastq.gz")
    message:
        "Linking file to output directory"
    shell:
        """
        ln -sr {input} {output}
        """

rule demux_bx:
    input:
        expand(outdir + "DATA_{IR}{ext}.fastq.gz", IR = ["R","I"], ext = [1,2])
    output:
        expand(outdir + inprefix + "_R{ext}.fastq.gz", ext = [1,2])
    message:
        "Moving barcodes into read headers"
    params:
        prfx = inprefix,
        outprfx = outdir + inprefix
    shell:
        """
        gen1.demux.o {params}
        """

rule split_samples_fw:
    input:
        outdir + inprefix + "_R1.fastq.gz"
    output:
        outdir + "{sample}.F.fq.gz"
    message:
        "Demultiplexing forward reads: {wildcards.sample}"
    params:
        parameters = paramspace.instance
    shell:
        """
        i="{wildcards.barcode}"
        zgrep A.."$i"B..D -A3 {input} | grep -v "^\-\-$"  | gzip > {output}
        """

rule split_samples_rv:
    input:
        outdir + inprefix + "_R2.fastq.gz"
    output:
        outdir + "{sample}.F.fq.gz"
    message:
        "Demultiplexing reverse reads: {wildcards.sample}"
    params:
        parameters = paramspace.instance
    shell:
        """
        i="{wildcards.barcode}"
        zgrep A.."$i"B..D -A3 {input} | grep -v "^\-\-$"  | gzip > {output}
        """

rule all
    input:
expand("Impute/{stitchparams}/variants.imputed.bcf", stitchparams=paramspace.instance_patterns),

        fw = ,
        rv = 
    default_target: True
    message:
        "Demultiplexing has finished!"    




        """
        for i in {01..96}; do
            zgrep A..C"$i"B..D -A3 | grep -v "^\-\-$"  > {output}
        done
        """