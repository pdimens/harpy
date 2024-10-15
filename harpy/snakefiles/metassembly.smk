containerized: "docker://pdimens/harpy:latest"

import os
import logging

onstart:
    logger.logger.addHandler(logging.FileHandler(config["snakemake_log"]))
onsuccess:
    os.remove(logger.logfile)
onerror:
    os.remove(logger.logfile)

FQ1 = config["inputs"]["fastq_r1"],
FQ2 = config["inputs"]["fastq_r2"]
outdir = config["output_directory"]
envdir = os.getcwd() + "/.harpy_envs"
max_mem = config["spades"]["max_memory"]
k_param = config["spades"]["k"]
extra = config["spades"].get("extra", "") 

rule sort_by_barcode:
    input:
        fq_f = FQ1,
        fq_r = FQ2
    output:
        fq_f = temp(f"{outdir}/fastq_preproc/tmp.R1.fq"),
        fq_r = temp(f"{outdir}/fastq_preproc/tmp.R2.fq")
    params:
        barcode_tag = config["barcode_tag"].upper()
    container:
        None
    shell:
        """
        samtools import -T "*" {input} |
        samtools sort -O SAM -t {params.barcode_tag} |
        samtools fastq -T "*" -1 {output.fq_f} -2 {output.fq_r}
        """
        
rule format_barcode:
    input:
        f"{outdir}/fastq_preproc/tmp.R{{FR}}.fq"
    output:
        temp(f"{outdir}/fastq_preproc/input.R{{FR}}.fq.gz")
    params:
        config["barcode_tag"].upper()
    container:
        None
    shell:
        "sed 's/{params}:Z:[^[:space:]]*/&-1/g' {input} | bgzip > {output}"

rule error_correction:
    input:
        FQ_R1 = outdir + "/fastq_preproc/input.R1.fq.gz",
        FQ_R2 = outdir + "/fastq_preproc/input.R2.fq.gz"
    output:
        outdir + "/error_correction/corrected/input.R1.fq00.0_0.cor.fastq.gz",
        outdir + "/error_correction/corrected/input.R2.fq00.0_0.cor.fastq.gz",
        outdir + "/error_correction/corrected/input.R_unpaired00.0_0.cor.fastq.gz"
    params:
        outdir = outdir + "/error_correction",
        k = k_param,
        mem = max_mem // 1000,
        extra = extra
    log:
        outdir + "/logs/error_correct.log"
    conda:
        f"{envdir}/assembly.yaml"
    threads:
        workflow.cores
    resources:
        mem_mb=max_mem
    shell:
        "metaspades.py -t {threads} -m {params.mem} -k {params.k} {params.extra} -1 {input.FQ_R1} -2 {input.FQ_R2} -o {params.outdir} --only-error-correction > {log}"

rule spades_assembly:
    input:
        FQ_R1C = outdir + "/error_correction/corrected/input.R1.fq00.0_0.cor.fastq.gz",
        FQ_R2C = outdir + "/error_correction/corrected/input.R2.fq00.0_0.cor.fastq.gz",
        FQ_UNC = outdir + "/error_correction/corrected/input.R_unpaired00.0_0.cor.fastq.gz"
    output:
        f"{outdir}/metaspades_assembly/contigs.fasta" 
    params:
        outdir = outdir + "/metaspades_assembly",
        k = k_param,
        mem = max_mem // 1000,
        extra = extra
    log:
        outdir + "/logs/metaspades_assembly.log"
    conda:
        f"{envdir}/assembly.yaml"
    threads:
        workflow.cores
    resources:
        mem_mb=max_mem
    shell:
        "metaspades.py -t {threads} -m {params.mem} -k {params.k} {params.extra} -1 {input.FQ_R1C} -2 {input.FQ_R2C} -s {input.FQ_UNC} -o {params.outdir} --only-assembler > {log}"

rule bwa_index:
    input:
        f"{outdir}/metaspades_assembly/contigs.fasta" 
    output:
        multiext(f"{outdir}/metaspades_assembly/contigs.fasta.", "ann", "bwt", "pac", "sa", "amb") 
    log:
        f"{outdir}/logs/bwa.index.log"
    conda:
        f"{envdir}/align.yaml"
    shell:
        "bwa index {input}"

rule bwa_align:
    input:
        multiext(f"{outdir}/metaspades_assembly/contigs.fasta.", "ann", "bwt", "pac", "sa", "amb"),
        fastq   = collect(outdir + "/fastq_preproc/input.R{X}.fq.gz", X = [1,2]),
        contigs = f"{outdir}/metaspades_assembly/contigs.fasta"
    output:
        f"{outdir}/reads-to-metaspades.bam"
    log:
        bwa = f"{outdir}/logs/align.bwa.log",
        samsort = f"{outdir}/logs/sort.alignments.log"
    threads:
        workflow.cores
    conda:
        f"{envdir}/align.yaml"
    shell:
        "bwa mem -C -t {threads} {input.contigs} {input.fastq} 2> {log.bwa} | samtools sort -O bam -o {output} - 2> {log.samsort}"

rule index_alignment:
    input:
        f"{outdir}/reads-to-metaspades.bam"
    output:
       f"{outdir}/reads-to-metaspades.bam.bai"
    log:
        f"{outdir}/logs/index.alignments.log"
    container:
        None
    shell:
        "samtools index {input} 2> {log}"

rule interleave_fastq:
    input:
        collect(outdir + "/fastq_preproc/input.R{FR}.fq.gz", FR = [1,2])
    output:
        f"{outdir}/fastq_preproc/interleaved.fq"
    container:
        None
    shell:
        "seqtk mergepe {input} > {output}"

rule athena_config:
    input:
        f"{outdir}/reads-to-metaspades.bam.bai",
        fastq = f"{outdir}/fastq_preproc/interleaved.fq",
        bam = f"{outdir}/reads-to-metaspades.bam",
        contigs = f"{outdir}/metaspades_assembly/contigs.fasta"
    output:
        f"{outdir}/athena/athena.config"
    params:
        threads = workflow.cores
    run:
        import json
        config_data = {
            "input_fqs": input.fastq,  
            "ctgfasta_path": input.contigs,  
            "reads_ctg_bam_path": input.bam,  
            "cluster_settings": {  
                "processes": params.threads,  
                "cluster_options": {  
                    "extra_params": {"run_local": "True"}  
                }  
            }  
        }  
        with open(output[0], "w") as conf:  
            json.dump(config_data, conf, indent=4)  

rule athena:
    input:
        multiext(f"{outdir}/reads-to-metaspades.", "bam", "bam.bai"),
        f"{outdir}/fastq_preproc/interleaved.fq",
        f"{outdir}/metaspades_assembly/contigs.fasta",
        config = f"{outdir}/athena/athena.config"
    output:
        temp(directory(collect(outdir + "/athena/{X}", X = ["results", "logs", "working"]))),
        f"{outdir}/athena/flye-input-contigs.fa",
        f"{outdir}/athena/athena.asm.fa",
    log:
        f"{outdir}/logs/athena.log"
    params:
        local_asm = f"{outdir}/athena/results/olc/flye-input-contigs.fa",
        final_asm = f"{outdir}/athena/results/olc/athena.asm.fa",
        result_dir = f"{outdir}/athena"
    conda:
        f"{envdir}/metassembly.yaml"
    shell:
        """
        athena-meta --config {input.config} 2> {log} &&\\
        mv {params.local_asm} {params.final_asm} {params.result_dir}      
        """

rule workflow_summary:
    default_target: True
    input:
        f"{outdir}/athena/athena.asm.fa"
    params:
        bx = config["barcode_tag"].upper(),
        k_param = k_param,
        max_mem = max_mem,
        extra = extra
    run:
        summary_template = f"""  
The harpy metassembly workflow ran using these parameters:  

FASTQ inputs were sorted by their linked-read barcodes:  
    samtools import -T "*" FQ1 FQ2 |  
    samtools sort -O SAM -t {{params.bx}} |  
    samtools fastq -T "*" -1 FQ_out1 -2 FQ_out2  

Barcoded-sorted FASTQ files had "-1" appended to the barcode to make them Athena-compliant:  
    sed 's/{{params.bx}}:Z:[^[:space:]]*/&-1/g' FASTQ | bgzip > FASTQ_OUT  

Reads were assembled using metaspades:  
    k values: {{params.k_param}}  
    maximum memory: {{params.max_mem}}  
    extra parameters: {{params.extra}}  

Original input FASTQ files were aligned to the metagenome using BWA:  
    bwa mem -C -p metaspades.contigs FQ1 FQ2 | samtools sort -O bam -  

Barcode-sorted Athena-compliant sequences were interleaved with seqtk:  
    seqtk mergepe FQ1 FQ2 > INTERLEAVED.FQ  

Athena ran with the config file Harpy built from the files created from the previous steps:  
    athena-meta --config athena.config  

The Snakemake workflow was called via command line:  
    {config["workflow_call"]}
"""
        with open(outdir + "/workflow/metassembly.summary", "w") as f:  
            f.write(summary_template.format(params=params))
