import os
import re
import gzip
import pysam
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
        temp(f"{outdir}/input.bam")
    threads:
        workflow.cores
    shell:
        """
        samtools import -@ {threads} -T "*" {input} -O BAM > {output}
        """

rule sample_barcodes:
    input:
        f"{outdir}/input.bam" if is_fastq else inputs[0]
    output:
        f"{outdir}/{prefix}.barcodes"
    log:
        f"{outdir}/logs/sampled_barcodes.log"
    threads:
        workflow.cores
    params:
        inv_prop = f"-i {invalids}",
        downsample_amt = f"-d {downsample}",
        bx_tag = "-b BX",
        random_seed = f"-r {random_seed}" if random_seed else ""
    shell:
        "extract_bxtags.py {params} {input} > {output} 2> {log}"

rule index_barcodes:
    input:
        lambda wc: infiles[wc.inputfile]
    output:
        bci = temp("{inputfile}.bci"),
        bai = temp("{inputfile}.bai") if not is_fastq else [],
        gzi = temp("{inputfile}i") if is_fastq and is_gzip else []
    threads:
        workflow.cores
    run:
        if is_fastq:
            gz_arg = "--gzip" if is_gzip else ""
            shell(f"LRez index fastq -t {threads} -f {input} -o {output.bci} {gz_arg}")
        else:
            shell(f"samtools index {input}")
            shell(f"LRez index bam --offsets -t {threads} -b {input} -o {output.bci}")

rule downsample:
    input:
        bam = inputs[0],
        bai = inputs[0] + ".bai",
        bc_index = inputs[0] + ".bci",
        bc_list = f"{outdir}/{prefix}.barcodes"
    output:
        sam = temp(f"{outdir}/{prefix}.sam"),
        bam = f"{outdir}/{prefix}.bam"
    threads:
        workflow.cores
    shell:
        """
        samtools view -H {input.bam} > {output.sam}
        LRez query bam -t {threads} -b {input.bam} -i {input.bc_index} -l {input.bc_list} >> {output.sam}
        samtools view -O BAM {output.sam} > {output.bam}
        """

rule downsample_read_1:
    input:
        fastq = inputs[0],
        bc_index = inputs[0] + ".bci",
        fq_index = inputs[0] + "i",
        bc_list = f"{outdir}/{prefix}.barcodes"
    output:
        f"{outdir}/{prefix}.R1.fq.gz"
    params:
        "--gzip" if is_gzip else ""
    threads:
        workflow.cores
    shell:
        "LRez query fastq -t {threads} -f {input.fastq} -i {input.bc_index} -l {input.bc_list} {params} | bgzip > {output}"

rule downsample_read_2:
    input:
        file = inputs[-1],
        bc_index = inputs[-1] + ".bci",
        fq_index = inputs[-1] + "i",
        bc_list = f"{outdir}/{prefix}.barcodes"
    output:
        f"{outdir}/{prefix}.R2.fq.gz"
    params:
        "--gzip" if is_gzip else ""
    threads:
        workflow.cores
    shell:
        "LRez query fastq -t {threads} -f {input.file} -i {input.bc_index} -l {input.bc_list} {params} | bgzip > {output}"

rule workflow_summary:
    default_target: True
    input:
        f"{outdir}/{prefix}.bam" if not is_fastq else [],
        collect(f"{outdir}/{prefix}.R" + "{FR}.fq.gz", FR = [1,2]) if is_fastq else []
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
        extraction += f"\textract_bxtags.py -i {invalids} -b BX -d {downsample} {params.random_seed} input.bam"
        summary.append(extraction)
        lrez = "The inputs were indexed and downsampled using LRez:\n"
        if is_fastq:
            lrez += "\tLRez query fastq -f fastq -i index.bci -l barcodes.txt"
        else:
            lrez += "\tLRez query bam -b bam -i index.bci -l barcodes.txt"
        summary.append(lrez)
        with open(outdir + "/workflow/downsample.summary", "w") as f:
            f.write("\n\n".join(summary))
