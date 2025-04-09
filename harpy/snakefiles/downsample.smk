import os
import re
import gzip
import pysam
import logging

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake_log"])
    logger.addHandler(logfile_handler)

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
        temp("input.bam")
    threads:
        workflow.cores
    shell:
        """
        samtools import -@ {threads} -T "*" {input} -O BAM > {output}
        """

rule sample_barcodes:
    input:
        "input.bam" if is_fastq else inputs[0]
    output:
        f"{prefix}.barcodes"
    log:
        "logs/sampled_barcodes.log"
    threads:
        workflow.cores
    params:
        inv_prop = f"-i {invalids}",
        downsample_amt = f"-d {downsample}",
        bx_tag = "-b BX",
        random_seed = f"-r {random_seed}" if random_seed else ""
    shell:
        "extract_bxtags.py {params} {input} > {output} 2> {log}"

if is_fastq:
    rule index_input:
        input:
            lambda wc: infiles[wc.inputfile]
        output:
            bci = temp("{inputfile}.bci"),
            gzi = temp("{inputfile}i") if is_fastq and is_gzip else []
        log:
            "logs/{inputfile}.index.log"
        params:
            gz_arg = "--gzip" if is_gzip else ""
        threads:
            workflow.cores
        shell:
            "LRez index fastq -t {threads} -f {input} -o {output.bci} {params.gz_arg}"
else:
    rule index_input:
        input:
            inputs[0]
        output:
            inputs[0] + ".bai"
        log:
            "log/" + os.path.basename(inputs[0]) + ".index.log" 
        shell:
            "samtools index {input}"

rule downsample:
    input:
        bam = inputs[0],
        bai = inputs[0] + ".bai",
        bc_list = f"{prefix}.barcodes"
    output:
        bam = f"{prefix}.bam"
    threads:
        workflow.cores
    shell:
        "samtools view -O BAM -h -D BX:{input.bc_list} {input.bam} > {output.bam}"

rule downsample_read_1:
    input:
        fastq = inputs[0],
        bc_index = inputs[0] + ".bci",
        fq_index = inputs[0] + "i",
        bc_list = f"{prefix}.barcodes"
    output:
        f"{prefix}.R1.fq.gz"
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
        bc_list = f"{prefix}.barcodes"
    output:
        f"{prefix}.R2.fq.gz"
    params:
        "--gzip" if is_gzip else ""
    threads:
        workflow.cores
    shell:
        "LRez query fastq -t {threads} -f {input.file} -i {input.bc_index} -l {input.bc_list} {params} | bgzip > {output}"

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
        extraction += f"\textract_bxtags.py -i {invalids} -b BX -d {downsample} {params.random_seed} input.bam"
        summary.append(extraction)
        lrez = "The inputs were indexed and downsampled using:\n"
        if is_fastq:
            lrez += "\tLRez query fastq -f fastq -i index.bci -l barcodes.txt"
        else:
            lrez += "\tsamtools view -O BAM -h -D BX:barcodes.txt input.bam"
        summary.append(lrez)
        with open("workflow/downsample.summary", "w") as f:
            f.write("\n\n".join(summary))
