import os
import re
import sys
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

outdir        = config["output_directory"]
envdir        = os.path.join(os.getcwd(), outdir, "workflow", "envs")
inputs        = config["inputs"]
invalids      = config["invalid_strategy"]
infiles       = dict(zip(inputs, inputs))
is_fastq      = True if len(inputs) == 2 else False
random_seed = config.get("random_seed", None)
rng = random.Random(random_seed) if random_seed else random.Random()

if filetype == "fq":
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
    threads:
        workflow.cores
    params:
        inv_prop = invalid_proportion,
        downsample_amt = downsample
    run:
        invalid_pattern = re.compile(r'[AaBbCcDd]00')
        barcodes = set()
        with pysam.AlignmentFile(sorted_bam, "rb", check_sq=False) as infile"
            for record in infile:
                try:
                    barcode = record.get_tag(bx_tag)
                    if invalid_pattern.search(barcode):
                        # invalid barcode retention
                        if rng.random() > params.inv_prop:
                            continue
                except KeyError:
                    continue
                barcodes.add(barcode)
            n_bc = len(barcodes)
            if n_bc < params.downsample_amt:
                raise ValueError(f"The input has fewer barcodes ({n_bc}) than the requested downsampling amount ({params.downsample_amt})")
            else:
                with open(output[0], "w") as bc_out:
                    _= [bc_out.write(f"{i}\n" for i in rng.sample(sorted(barcodes), params.downsample_amt))]

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
        f"{outdir}/downsample.bam"
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
        f"{outdir}/downsample.R1.fq.gz"
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
        f"{outdir}/downsample.R2.fq.gz"
    params:
        "--gzip" if is_gzip else ""
    threads:
        workflow.cores
    shell:
        "LRez query fastq -t {threads} -f {input.file} -i {input.bc_index} -l {input.bc_list} {params} > {output}"


rule workflow_summary:
    input:
        f"{outdir}/downsample.bam" if not is_fastq else [],
        collect(f"{outdir}/downsample.R" + "{FR}.fq.gz", FR = [1,2]) if is_fastq else []
