import os
import re
import sys
import logging as pylogging

rule metaspades:
    input:
    output:
        F_fq = f"{outidr}/metaspades/corrected/",
        R_fq = f"{outidr}/metaspades/corrected/",
        spades_contigs = f"{outidr}/metaspades/contigs.fasta" 
    log:
        f"{outdir}/logs/metaspades.log"
    params:
        f"{outdir}/metaspades"
    shell:
        "metaspades.py --12 {input} -o {params} 2> {log}"

rule bwa_index:
    input:
        f"{outidr}/metaspades/contigs.fasta" 
    output:
        multiext(f"{outidr}/metaspades/contigs.fasta.", "ann", "bwt", "pac", "sa", "amb") 
    log:
        f"{outdir}/logs/bwa.index.log"
    shell:
        "bwa index {input}"

rule bwa_align:
    input:
        fastq   = path/to/fastq,
        contigs = f"{outidr}/metaspades/contigs.fasta",
        indices = multiext(f"{outidr}/metaspades/contigs.fasta.", "ann", "bwt", "pac", "sa", "amb")
    output:
        f"{outdir}/align/reads-to-metaspades.bam"
    log:
        bwa = f"{outdir}/logs/align.bwa.log",
        samsort = f"{outdir}/logs/sort.alignments.log"
    shell:
        "bwa mem -C -p {input.contigs} /path/to/reads 2> {log.bwa} | samtools sort -O bam -o {output} - 2> {log.samsort}"

rule index_alignment:
    input:
        f"{outdir}/align/reads-to-metaspades.bam"
    output:
`       f"{outdir}/align/reads-to-metaspades.bam.bai"
    log:
        f"{outdir}/logs/index.alignments.log"
    shell:
        "samtools index {input} 2> {log}"

rule athena_config:
    input:
        fastq = path/to/fasq
        bam = f"{outdir}/align/reads-to-metaspades.bam",
        bai = f"{outdir}/align/reads-to-metaspades.bam.bai",
        contigs = f"{outidr}/metaspades/contigs.fasta"
    output:
        f"{outdir}/athena/athena.config"
    run:
        with open(output[0], "w") as conf:
        _ = conf.write("{\n")
        _ = conf.write("\"input_fqs\": \"/path/to/fq\",\n")
        _ = conf.write("\"ctgfasta_path\": \"/path/to/seeds.fa\",\n")
        _ = conf.write("\"reads_ctg_bam_path\": \"/path/to/reads_2_seeds.bam\"\n")
        _ = conf.write("}\n")

rule athena:
    input:
        fastq = path/to/fasq,
        bam = f"{outdir}/align/reads-to-metaspades.bam",
        bai = f"{outdir}/align/reads-to-metaspades.bam.bai",
        contigs = f"{outidr}/metaspades/contigs.fasta",
        config = f"{outdir}/athena/athena.config"
    output:
        local_asm = f"{outdir}/athena/results/olc/flye-input-contigs.fa",
        final_asm = f"{outdir}/athena/results/olc/athena.asm.fa"
    log:
        f"{outdir}/logs/athena.log"
    shell:

rule pangaea:
    input:
        F_fq = f"{outidr}/metaspades/corrected/",
        R_fq = f"{outidr}/metaspades/corrected/",
        spades_contigs = f"{outidr}/metaspades/corrected/contigs.fasta" 
    output:

    log:
        f"{outdir}/logs/pangaea.log"
    shell:
        "python pangaea_path/pangaea.py -1 {input.F_fq} -2 {input.R_fq} -sp {input.spades_contigs} -lc {athena_local} -at {athena_hybrid} -lt 10,30 -c 30 -o {outdir} 2> {log}"

rule workflow_summary:
