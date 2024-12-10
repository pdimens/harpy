import os
import re
import gzip
import pysam
import random
import logging

onstart:
    logger.logger.addHandler(logging.FileHandler(config["snakemake_log"]))
onsuccess:
    os.remove(logger.logfile)
onerror:
    os.remove(logger.logfile)

outdir      = config["output_directory"]
inputs      = config["inputs"]
invalids    = config["invalid_proportion"]
random_seed = config.get("random_seed", None)
downsample  = config["downsample"]
prefix      = config["prefix"]
infiles     = dict(zip(inputs, inputs))
is_fastq    = True if len(inputs) == 2 else False
rng         = random.Random(random_seed) if random_seed else random.Random()

if is_fastq:
    # determine if a file is gzipped
    try:
        with gzip.open(inputs[0], 'rt') as f:
            f.read(10)
        is_gzip = True
    except gzip.BadGzipFile:
        is_gzip = False
else:
    is_gzip = False

rule bam_convert:
    input:
        inputs
    output:
        temp(f"{outdir}/input.bam")
    threads:
        workflow.cores
    shell:
        """
        samtools import -@ {threads} -T "*" {input}
        """

rule extract_barcodes:
    input:
        f"{outdir}/input.bam" if is_fastq else inputs[0]
    output:
        f"{outdir}/sampled_barcodes.txt"
    log:
        f"{outdir}/logs/sampled_barcodes.log"
    threads:
        workflow.cores
    params:
        inv_prop = f"-i {invalid_proportion}",
        downsample_amt = f"-d {downsample}",
        bx_tag = "-b BX"
    shell:
        "extract_bxtags.py {params} {input} > {output} 2> {log}"

rule index_barcodes:
    input:
        lambda wc: infiles[wc.inputfile]
    output:
        "{inputfile}.bci"
    threads:
        workflow.cores
    run:
        if is_fastq:
            gz_arg = "--gzip" if is_gzip else ""
            shell(f"LRez index fastq -t {threads} -f {input} -o {output} {gz_arg}")
        else:
            shell(f"LRez index bam --offsets -t {threads} -b {input} -o {output}")

rule downsample:
    input:
        file = infiles[0],
        bc_index = infiles[0] + ".bci",
        bc_list = f"{outdir}/sampled_barcodes.txt"
    output:
        f"{outdir}/{prefix}.bam"
    threads:
        workflow.cores
    shell:
        "LRez query bam -t {threads} -bam {input.file} -i {input.bc_index} -l {input.bc_list} > {output}"

rule downsample_read_1:
    input:
        file = inputs[0],
        bc_index = inputs[0] + ".bci",
        bc_list = f"{outdir}/sampled_barcodes.txt"
    output:
        f"{outdir}/{prefix}.R1.fq.gz"
    params:
        "--gzip" if is_gzip else ""
    threads:
        workflow.cores
    shell:
        "LRez query fastq -t {threads} -f {input.file} -i {input.bc_index} -l {input.bc_list} {params} > {output}"
 
rule downsample_read_2:
    input:
        file = inputs[1],
        bc_index = inputs[1] + ".bci",
        bc_list = f"{outdir}/sampled_barcodes.txt"
    output:
        f"{outdir}/{prefix}.R2.fq.gz"
    params:
        "--gzip" if is_gzip else ""
    threads:
        workflow.cores
    shell:
        "LRez query fastq -t {threads} -f {input.file} -i {input.bc_index} -l {input.bc_list} {params} > {output}"

rule workflow_summary:
    input:
        f"{outdir}/{prefix}.bam" if not is_fastq else [],
        collect(f"{outdir}/{prefix}.R" + "{FR}.fq.gz", FR = [1,2]) if is_fastq else []
    run:
        summary = ["The harpy downsample workflow ran using these parameters:"]
        summary.append(f"The provided input file(s):\n" + "\n\t".join(inputs))
        convs = "The FASTQ files were converted into a BAM file with:\n"
        convs += "\tsamtools import -T * fastq1 fastq2"
        if is_fastq:
            summary.append(convs)
        extraction = "Barcodes were extracted and sampled using:\n"
        extraction += f"\textract_bxtags.py -i {invalids} -b BX -d {downsample} input.bam"
        summary.append(extraction)
        lrez = "The inputs were indexed and downsampled using LRez:\n"
        if is_fastq:
            lrez += "\tLRez query fastq -f fastq -i index.bci -l barcodes.txt"
        else:
            lrez += "\tLRez query bam -b bam -i index.bci -l barcodes.txt"
        summary.append(lrez)
        with open(outdir + "/workflow/downsample.summary", "w") as f:
            f.write("\n\n".join(summary))
