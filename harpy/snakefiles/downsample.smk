import os
import re
import gzip
import pysam
import logging

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake"]["log"])
    logger.addHandler(logfile_handler)

inputs      = config["inputs"]
invalids    = config["invalid_proportion"]
random_seed = config.get("random_seed", None)
downsample  = config["downsample"]
bc_tag      = config["barcode-tag"]
prefix      = config["prefix"]
infiles     = dict(zip(inputs, inputs))
is_fastq    = True if len(inputs) == 2 else False
bam_file = "input.bam" if is_fastq else inputs[0]

if is_fastq:
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
        temp("input.bam")
    threads:
        workflow.cores
    shell:
        """
        samtools import -@ {threads} -T "*" {input} -O BAM > {output}
        """

rule sample_barcodes:
    input:
        bam_file
    output:
        f"{prefix}.barcodes"
    log:
        "logs/sampled_barcodes.log"
    threads:
        workflow.cores
    params:
        inv_prop = f"-i {invalids}",
        downsample_amt = f"-d {downsample}",
        bx_tag = f"-b {bc_tag}",
        random_seed = f"-r {random_seed}" if random_seed else ""
    shell:
        "extract_bxtags {params} {input} > {output} 2> {log}"

rule index_input:
    input:
        bam_file
    output:
        temp(bam_file + ".bai")
    log:
        "logs/" + os.path.basename(bam_file) + ".index.log" 
    shell:
        "samtools index {input}"

rule downsample:
    input:
        bam = bam_file,
        bai = bam_file + ".bai",
        bc_list = f"{prefix}.barcodes"
    output:
        bam = temp(f"{prefix}.bam") if is_fastq else f"{prefix}.bam"
    params:
        bc_tag
    threads:
        workflow.cores
    shell:
        "samtools view -O BAM -h -D {params}:{input.bc_list} {input.bam} > {output.bam}"

rule revert_to_fastq:
    input:
        f"{prefix}.bam"
    output:
        R1 = f"{prefix}.R1.fq.gz",
        R2 = f"{prefix}.R2.fq.gz"
    threads:
        workflow.cores
    shell:
        "samtools fastq -@ {threads} -T \"*\" -1 {output.R1} -2 {output.R2} {input}"

rule workflow_summary:
    default_target: True
    input:
        f"{prefix}.bam" if not is_fastq else [],
        collect(f"{prefix}.R" + "{FR}.fq.gz", FR = [1,2]) if is_fastq else []
    params:
        random_seed = f"-r {random_seed}" if random_seed else ""
    run:
        summary = ["The harpy downsample workflow ran using these parameters:"]
        summary.append(f"The provided input file(s):\n\t" + "\n\t".join(inputs))
        convs = "The FASTQ files were converted into a BAM file with:\n"
        convs += "\tsamtools import -T * fastq1 fastq2"
        if is_fastq:
            summary.append(convs)
        extraction = "Barcodes were extracted and sampled using:\n"
        extraction += f"\textract_bxtags -i {invalids} -b BX -d {downsample} {params.random_seed} input.bam"
        summary.append(extraction)
        downsampled = "The inputs were indexed and downsampled using:\n"
        downsampled += f"\tsamtools view -O BAM -h -D {bc_tag}:barcodes.txt input.bam"
        summary.append(downsampled)
        revs = "The input fastq fles were reverted to FASTQ format with:\n"
        revs += "\tsamtools fastq -T \"*\" -1 R1 -2 R2 input.bam"
        if is_fastq:
            summary.append(revs)
        with open("workflow/downsample.summary", "w") as f:
            f.write("\n\n".join(summary))
