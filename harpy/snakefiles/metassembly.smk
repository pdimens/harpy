containerized: "docker://pdimens/harpy:latest"

import os
import logging

onstart:
    logger.logger.addHandler(logging.FileHandler(config["snakemake_log"]))
onsuccess:
    os.remove(logger.logfile)
onerror:
    os.remove(logger.logfile)

#cont_cov = config["contig_coverage"]
#clusters = config["clusters"]
outdir = config["output_directory"]
envdir = os.getcwd() + "/.harpy_envs"
max_mem = config["metaspades"]["max_memory"]
k_param = config["metaspades"]["k"]
extra = config["metaspades"].get("extra", "") 

rule sort_by_barcode:
    input:
        fq_f = config["inputs"]["fastq_r1"],
        fq_r = config["inputs"]["fastq_r2"]
    output:
        fq_f = temp(f"{outdir}/fastq_preproc/tmp.R1.fq"),
        fq_r = temp(f"{outdir}/fastq_preproc/tmp.R2.fq")
    params:
        config["barcode_tag"].upper()
    container:
        None
    shell:
        """
        samtools import -T "*" {input} |
        samtools sort -O SAM -t {params} |
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

rule metaspades:
    input:
        reads = collect(outdir + "/fastq_preproc/input.R{X}.fq.gz", X = [1,2])
    output:
        contigs = outdir + "/metaspades/contigs.fasta",
        scaffolds = outdir + "/metaspades/scaffolds.fasta",
        dir = directory(outdir + "/metaspades/intermediate_files"),
        corrected_F = outdir + "/metaspades/intermediate_files/corrected/input.R100.0_0.cor.fastq.gz",
        corrected_R = outdir + "/metaspades/intermediate_files/corrected/input.R200.0_0.cor.fastq.gz"
    params:
        k = k_param,
        extra = extra
    log:
        outdir + "/logs/spades.log"
    threads:
        workflow.cores
    resources:
        mem_mb=max_mem
    wrapper:
        "v4.7.1/bio/spades/metaspades"

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
        multiext(f"{outdir}/metaspades/contigs.fasta.", "ann", "bwt", "pac", "sa", "amb"),
        fastq   = collect(outdir + "/fastq_preproc/input.R{X}.fq.gz", X = [1,2]),
        contigs = f"{outdir}/metaspades/contigs.fasta"
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
        f"{outdir}/fastq_preproc/inverleaved.fq"
    container:
        None
    shell:
        "seqtk mergepe {input} > {output}"

rule athena_config:
    input:
        f"{outdir}/reads-to-metaspades.bam.bai",
        fastq = f"{outdir}/fastq_preproc/inverleaved.fq",
        bam = f"{outdir}/reads-to-metaspades.bam",
        contigs = f"{outdir}/metaspades/contigs.fasta"
    output:
        f"{outdir}/athena/athena.config"
    run:
        with open(output[0], "w") as conf:
            _ = conf.write("{\n")
            _ = conf.write(f"\"input_fqs\": \"{input.fastq}\",\n")
            _ = conf.write(f"\"ctgfasta_path\": \"{input.contigs}\",\n")
            _ = conf.write(f"\"reads_ctg_bam_path\": \"{input.bam}\"\n")
            _ = conf.write("}\n")

rule athena:
    input:
        multiext(f"{outdir}/reads-to-metaspades.", "bam", "bam.bai"),
        f"{outdir}/workflow/input.fq",
        f"{outdir}/metaspades/contigs.fasta",
        config = f"{outdir}/athena/athena.config"
    output:
        local_asm = f"{outdir}/athena/olc/flye-input-contigs.fa",
        final_asm = f"{outdir}/athena/olc/athena.asm.fa"
    log:
        f"{outdir}/logs/athena.log"
    conda:
        f"{envdir}/metassembly.yaml"
    shell:
        "athena-meta --config {input.config}"

#TODO figure this part out
# is it sorted reads, interleaved? Can I just use the starting ones? maybe just the corrected metaspades ones
#rule pangaea:
#    input:
#        F_fq = f"{outdir}/metaspades/corrected/corrected/input_100.0_0.cor.fastq.gz",
#        R_fq = f"{outdir}/metaspades/corrected/corrected/input_200.0_0.cor.fastq.gz",
#        spades_contigs = f"{outdir}/metaspades/corrected/contigs.fasta",
#        athena_local = f"{outdir}/athena/results/olc/flye-input-contigs.fa",
#        athena_hybrid = f"{outdir}/athena/results/olc/athena.asm.fa"
#    output:
#        f"{outdir}/pangaea/final.asm.fa"
#    params:
#        outdir = f"{outdir}/pangaea",
#        lt = cont_cov,
#        c =  clusters
#    log:
#        f"{outdir}/logs/pangaea.log"
#    shell:
#        "python pangaea_path/pangaea.py -1 {input.F_fq} -2 {input.R_fq} -sp {input.spades_contigs} -lc {input.athena_local} -at {input.athena_hybrid} -lt {params.lt} -c {params.c} -o {params.outdir} 2> {log}"

rule workflow_summary:
    default_target: True
    input:
        f"{outdir}/athena/olc/athena.asm.fa"
    params:
        bx = config["barcode_tag"].upper(),
        k_param = k_param,
        max_mem = max_mem,
        extra = extra
        #cont_cov = cont_cov,
        #clusters = clsuters,
    run:
        with open(outdir + "/workflow/metassembly.summary", "w") as f:
            _ = f.write("The harpy metassembly workflow ran using these parameters:\n\n")
            _ = f.write("FASTQ inputs were sorted by their linked-read barcodes:\n")
            _ = f.write("    samtools import -T \"*\" FQ1 FQ2 |\n")
            _ = f.write(f"    samtools sort -O SAM -t {params.bx} |\n")
            _ = f.write("    samtools fastq -T "*" -1 FQ_out1 -2 FQ_out2\n")
            _ = f.write("Barcoded-sorted FASTQ files had \"-1\" appended to the barcode to make them Athena-compliant:\n")
            _ = f.write(f"    sed 's/{params.bx}:Z:[^[:space:]]*/&-1/g' FASTQ | bgzip > FASTQ_OUT\n")
            _ = f.write(f"Reads were assembled using metaspades:\n")
            _ = f.write(f"    k values: {params.k_param}\n")
            _ = f.write(f"    maximum memory: {params.max_mem}\n")
            _ = f.write(f"    extra parameters: {params.extra}\n")
            _ = f.write("Original input FASTQ files were aligned to the metagenome using BWA:\n")
            _ = f.write("    bwa mem -C -p metaspades.contigs FQ1 FQ2 | samtools sort -O bam -\n")
            _ = f.write("Barcode-sorted Athena-compliant sequences were interleaved with seqtk:\n")
            _ = f.write("    seqtk mergepe FQ1 FQ2 > INTERLEAVED.FQ\n")
            _ = f.write("Athena ran with the config file Harpy built from the files created from the previous steps:\n")
            _ = f.write("    athena-meta --config athena.config\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")