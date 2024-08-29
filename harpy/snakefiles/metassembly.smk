import os
import logging as pylogging

FASTQ1 = config["inputs"]["fastq"]
FASTQ2 = config["inputs"].get("fastq", None)
cont_cov = config["contig_coverage"]
clusters = config["clusters"]
snakemake_log = config["snakemake_log"]
outdir = config["output_directory"]
envdir = os.getcwd() + "/.harpy_envs"

onstart:
    extra_logfile_handler = pylogging.FileHandler(snakemake_log)
    logger.logger.addHandler(extra_logfile_handler)

rule sort_fastq:
    input:
        fq = FASTQ1,
        fq2 = FASTQ2 if FASTQ2 else []
    output:
        temp(f"{outdir}/workflow/input.fq")
    params:
        config["barcode_tag"].upper()
    shell:
        """
        samtools import -T "*" {input} |
        samtools sort -O SAM -t {params} |
        samtools fastq -T "*" > {output}
        """

#TODO METASPADES NEEDS MORE PARAMETERS
rule metaspades:
    input:
        f"{outdir}/workflow/input.fq"
    output:
        F_fq = f"{outdir}/metaspades/corrected/corrected/input_1.00.0_0.cor.fastq.gz",
        R_fq = f"{outdir}/metaspades/corrected/corrected/input_2.00.0_0.cor.fastq.gz",
        spades_contigs = f"{outdir}/metaspades/contigs.fasta" 
    log:
        f"{outdir}/logs/metaspades.log"
    params:
        outdir = f"{outdir}/metaspades",
        continue_arg = "--restart-from last" if os.path.exists(f"{outdir}/metaspades/pipeline_state/") else "",
        mem_gb = lambda wildcards, resources: resources.mem_mb//1000,
        infile = "" if os.path.exists(f"{outdir}/metaspades/pipeline_state/") else f"--12 {outdir}/workflow/input.fq"
    threads:
        workflow.cores
    resources:
        mem_mb = 20000
    conda:
        f"{envdir}/metassembly.yaml"
    shell:
        "metaspades.py {params.infile} {params.continue_arg} -o {params.outdir} -t {threads} -m {params.mem_gb} > {log}"

rule all:
    default_target: True
    input:
        f"{outdir}/metaspades/contigs.fasta"

rule bwa_index:
    input:
        f"{outdir}/metaspades/contigs.fasta" 
    output:
        multiext(f"{outdir}/metaspades/contigs.fasta.", "ann", "bwt", "pac", "sa", "amb") 
    log:
        f"{outdir}/logs/bwa.index.log"
    conda:
        f"{envdir}/align.yaml"
    shell:
        "bwa index {input}"

rule bwa_align:
    input:
        fastq   = f"{outdir}/workflow/input.fq",
        contigs = f"{outdir}/metaspades/contigs.fasta",
        indices = multiext(f"{outdir}/metaspades/contigs.fasta.", "ann", "bwt", "pac", "sa", "amb")
    output:
        f"{outdir}/align/reads-to-metaspades.bam"
    log:
        bwa = f"{outdir}/logs/align.bwa.log",
        samsort = f"{outdir}/logs/sort.alignments.log"
    conda:
        f"{envdir}/align.yaml"
    shell:
        "bwa mem -C -p {input.contigs} /path/to/reads 2> {log.bwa} | samtools sort -O bam -o {output} - 2> {log.samsort}"

rule index_alignment:
    input:
        f"{outdir}/align/reads-to-metaspades.bam"
    output:
       f"{outdir}/align/reads-to-metaspades.bam.bai"
    log:
        f"{outdir}/logs/index.alignments.log"
    conda:
        f"{envdir}/align.yaml"
    shell:
        "samtools index {input} 2> {log}"

rule athena_config:
    input:
        fastq = f"{outdir}/workflow/input.fq",
        bam = f"{outdir}/align/reads-to-metaspades.bam",
        bai = f"{outdir}/align/reads-to-metaspades.bam.bai",
        contigs = f"{outdir}/metaspades/contigs.fasta"
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
        fastq = f"{outdir}/workflow/input.fq",
        bam = f"{outdir}/align/reads-to-metaspades.bam",
        bai = f"{outdir}/align/reads-to-metaspades.bam.bai",
        contigs = f"{outdir}/metaspades/contigs.fasta",
        config = f"{outdir}/athena/athena.config"
    output:
        local_asm = f"{outdir}/athena/results/olc/flye-input-contigs.fa",
        final_asm = f"{outdir}/athena/results/olc/athena.asm.fa"
    log:
        f"{outdir}/logs/athena.log"
    conda:
        f"{envdir}/metassembly.yaml"
    shell:
        "athena-meta --config {input.config}"

#TODO figure this part out
# is it sorted reads, interleaved? Can I just use the starting ones? maybe just the corrected metaspades ones
rule pangaea:
    input:
        F_fq = f"{outdir}/metaspades/corrected/corrected/input_1.00.0_0.cor.fastq.gz",
        R_fq = f"{outdir}/metaspades/corrected/corrected/input_2.00.0_0.cor.fastq.gz",
        spades_contigs = f"{outdir}/metaspades/corrected/contigs.fasta",
        athena_local = f"{outdir}/athena/results/olc/flye-input-contigs.fa",
        athena_hybrid = f"{outdir}/athena/results/olc/athena.asm.fa"
    output:
        f"{outdir}/pangaea/final.asm.fa"
    params:
        outdir = f"{outdir}/pangaea",
        lt = cont_cov,
        c =  clusters
    log:
        f"{outdir}/logs/pangaea.log"
    shell:
        "python pangaea_path/pangaea.py -1 {input.F_fq} -2 {input.R_fq} -sp {input.spades_contigs} -lc {input.athena_local} -at {input.athena_hybrid} -lt {params.lt} -c {params.c} -o {params.outdir} 2> {log}"

#rule workflow_summary:
