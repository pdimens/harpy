import os
import logging

onstart:
    logger.logger.addHandler(logging.FileHandler(config["snakemake_log"]))
onsuccess:
    os.remove(logger.logfile)
onerror:
    os.remove(logger.logfile)

FASTQ1 = config["inputs"]["fastq"]
FASTQ2 = config["inputs"].get("fastq2")
cont_cov = config["contig_coverage"]
clusters = config["clusters"]
outdir = config["output_directory"]
envdir = os.getcwd() + "/.harpy_envs"

rule barcode_sort:
    input:
        fq_f = FASTQ1,
        fq_r = FASTQ2
    output:
        fq_f = temp(f"{outdir}/workflow/input.R1.fq"),
        fq_r = temp(f"{outdir}/workflow/input.R2.fq")
    params:
        config["barcode_tag"].upper()
    shell:
        """
        samtools import -T "*" {input} |
        samtools sort -O SAM -t {params} |
        samtools fastq -T "*" -1 {output.fq_f} -2 {output.fq_r}
        sed -i 's/{params}:Z[^[:space:]]*/&-1/g' {output.fq_f}
        sed -i 's/{params}:Z[^[:space:]]*/&-1/g' {output.fq_r}
        """

rule metaspades:
    input:
        reads = collect(outdir + "/workflow/input.R{X}.fq", X = [1,2])
    output:
        contigs = outdir + "/metaspades/contigs.fasta",
        scaffolds = outdir + "/metaspades/scaffolds.fasta",
        dir = directory(outdir + "/metaspades/intermediate_files"),
        corrected_F = outdir + "/metaspades/intermediate_files/corrected/input.R100.0_0.cor.fastq.gz",
        corrected_R = outdir + "/metaspades/intermediate_files/corrected/input.R200.0_0.cor.fastq.gz"
    params:
        k="auto"
        #extra="--only-assembler",
    log:
        outdir + "/logs/spades.log",
    threads:
        workflow.cores
    resources:
        mem_mem=250000,
        time=60 * 24,
    wrapper:
        "v4.3.0/bio/spades/metaspades"

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
        fastq   = collect(outdir + "/workflow/input.R{X}.fq", X = [1,2]),
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
        "bwa mem -C -p {input.contigs} {input.fastq} 2> {log.bwa} | samtools sort -O bam -o {output} - 2> {log.samsort}"

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
    conda:
        f"{envdir}/metassembly.yaml"
    run:
        with open(output[0], "w") as conf:
            _ = conf.write("{\n")
            _ = conf.write("\"input_fqs\": \"/path/to/fq\",\n")
            _ = conf.write("\"ctgfasta_path\": \"/path/to/seeds.fa\",\n")
            _ = conf.write("\"reads_ctg_bam_path\": \"/path/to/reads_2_seeds.bam\"\n")
            _ = conf.write("}\n")

rule all:
    default_target: True
    input:
        f"{outdir}/athena/athena.config"

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
        F_fq = f"{outdir}/metaspades/corrected/corrected/input_100.0_0.cor.fastq.gz",
        R_fq = f"{outdir}/metaspades/corrected/corrected/input_200.0_0.cor.fastq.gz",
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
