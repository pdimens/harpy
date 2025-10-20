containerized: "docker://pdimens/harpy:latest"

import os
import logging

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake"]["log"])
    logger.addHandler(logfile_handler)

FQ1 = config["inputs"]["fastq_r1"]
FQ2 = config["inputs"]["fastq_r2"]
BX_TAG = config["linkedreads"]["barcode_tag"]
max_mem = config["spades"]["max_memory"]
k_param = config["spades"]["k"]
ignore_bx = config["spades"]["ignore_barcodes"]
extra = config["spades"].get("extra", "")
spadesdir = f"{'cloudspades' if not ignore_bx else 'spades'}_assembly"
force_athena = config["athena"]["force"]
skip_reports  = config["reports"]["skip"]
organism = config["reports"]["organism_type"]
lineage_map = {
    "eukaryote": "eukaryota",
    "fungus": "fungi",
    "bacteria": "bacteria"
}
lineagedb = lineage_map.get(organism, "bacteria")
odb_version = 12

rule sort_by_barcode:
    input:
        fq_f = FQ1,
        fq_r = FQ2
    output:
        fq_f = temp("fastq_preproc/tmp.R1.fq"),
        fq_r = temp("fastq_preproc/tmp.R2.fq")
    log:
        "logs/sort_by_barcode.log"
    params:
        barcode_tag = BX_TAG
    threads:
        workflow.cores
    shell:
        """
        {{
            samtools import -T "*" {input} |
            samtools sort -@ {threads} -O SAM -t {params.barcode_tag} |
            samtools fastq -T "*" -1 {output.fq_f} -2 {output.fq_r}
        }} 2> {log}
        """

rule format_barcode:
    input:
        "fastq_preproc/tmp.R{FR}.fq"
    output:
        temp("fastq_preproc/input.R{FR}.fq.gz")
    params:
        barcode_tag = BX_TAG
    shell:
        "sed 's/{params}:Z:[^[:space:]]*/&-1/g' {input} | bgzip > {output}"

rule error_correction:
    input:
        FQ_R1 = "fastq_preproc/input.R1.fq.gz",
        FQ_R2 = "fastq_preproc/input.R2.fq.gz"
    output:
        "error_correction/corrected/input.R1.fq00.0_0.cor.fastq.gz",
        "error_correction/corrected/input.R2.fq00.0_0.cor.fastq.gz",
        "error_correction/corrected/input.R_unpaired00.0_0.cor.fastq.gz"
    params:
        "--only-error-correction",
        "-o error_correction",
        f"-k {k_param}",
        f"-m {max_mem // 1000}",
        extra
    log:
        "logs/error_correct.log"
    threads:
        workflow.cores
    resources:
        mem_mb=max_mem
    conda:
        "envs/assembly.yaml"
    container:
        "docker://pdimens/harpy:assembly_latest"
    shell:
        "metaspades -t {threads} {params} -1 {input.FQ_R1} -2 {input.FQ_R2} > {log}"

rule spades_assembly:
    input:
        fastq_R1C = "error_correction/corrected/input.R1.fq00.0_0.cor.fastq.gz",
        fastq_R2C = "error_correction/corrected/input.R2.fq00.0_0.cor.fastq.gz",
        fastq_UNC = "error_correction/corrected/input.R_unpaired00.0_0.cor.fastq.gz"
    output:
        "spades_assembly/contigs.fasta" 
    params:
        "--only-assembler",
        "-o spades_assembly",
        f"-k {k_param}",
        f"-m {max_mem // 1000}",
        extra
    log:
        "logs/spades_assembly.log"
    threads:
        workflow.cores
    resources:
        mem_mb=max_mem
    conda:
        "envs/spades.yaml"
    container:
        None
    shell:
        "metaspades.py -t {threads} {params} -1 {input.fastq_R1C} -2 {input.fastq_R2C} -s {input.fastq_UNC} > {log}"

rule cloudspades_metassembly:
    input:
        fastq_R1 = FQ1,
        fastq_R2 = FQ2
    output:
        "cloudspades_assembly/contigs.fasta",
        "cloudspades_assembly/scaffolds.fasta"
    params:
        "--meta",
        f"-o {spadesdir}",
        f"-k {k_param}",
        f"-m {max_mem // 1000}",
        extra
    log:
        "logs/assembly.log"
    conda:
        "envs/assembly.yaml"
    threads:
        workflow.cores
    resources:
        mem_mb = max_mem
    shell:
        "spades.py -t {threads} {params} --gemcode1-1 {input.fastq_R1} --gemcode1-2 {input.fastq_R2} > {log}"

rule index_contigs:
    input:
        f"{spadesdir}/contigs.fasta"
    output:
        multiext(f"{spadesdir}/contigs.fasta.", "ann", "bwt", "pac", "sa", "amb") 
    log:
        "logs/bwa.index.log"
    conda:
        "envs/align.yaml"
    shell:
        "bwa index {input}"

rule align_to_contigs:
    input:
        multiext(f"{spadesdir}/contigs.fasta.", "ann", "bwt", "pac", "sa", "amb"),
        fastq   = collect("fastq_preproc/input.R{X}.fq.gz", X = [1,2]),
        contigs = f"{spadesdir}/contigs.fasta"
    output:
        bam = temp("reads-to-spades.bam"),
        bai = temp("reads-to-spades.bam.bai")
    log:
        "logs/align.to.contigs.log"
    params:
        f"-C -v 2 -t {workflow.cores - 1}"
    threads:
        workflow.cores
    conda:
        "envs/align.yaml"
    shell:
        """
        {{
            bwa mem {params} {input.contigs} {input.fastq} |
                samtools view -h -F 4 -q 5 |
                samtools sort -@ 1 -O bam --write-index -o {output.bam}##idx##{output.bai} -
        }} 2> {log}
        """

rule interleave_fastq:
    input:
        collect("fastq_preproc/input.R{FR}.fq.gz", FR = [1,2])
    output:
        temp("fastq_preproc/interleaved.fq")
    container:
        None
    shell:
        "seqtk mergepe {input} > {output}"

rule athena_config:
    input:
        "reads-to-spades.bam.bai",
        fastq = "fastq_preproc/interleaved.fq",
        bam = "reads-to-spades.bam",
        contigs = f"{spadesdir}/contigs.fasta"
    output:
        "athena/athena.config"
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

rule athena_metassembly:
    input:
        multiext("reads-to-spades.", "bam", "bam.bai"),
        "fastq_preproc/interleaved.fq",
        f"{spadesdir}/contigs.fasta",
        config = "athena/athena.config"
    output:
        temp(directory(collect("athena/{X}", X = ["results", "logs", "working"]))),
        "athena/flye-input-contigs.fa",
        "athena/athena.asm.fa",
    log:
        "logs/athena.log"
    params:
        force = "--force_reads" if force_athena else "",
        local_asm = "athena/results/olc/flye-input-contigs.fa",
        final_asm = "athena/results/olc/athena.asm.fa",
        result_dir = "athena"
    conda:
        "envs/metassembly.yaml"
    shell:
        """
        athena-meta {params.force} --config {input.config} &> {log} &&\\
        mv {params.local_asm} {params.final_asm} {params.result_dir}      
        """

rule QUAST_assessment:
    input:
        contigs = f"{spadesdir}/contigs.fasta",
        scaffolds = "athena/athena.asm.fa",
        fastq_f = FQ1,
        fastq_r = FQ2
    output:
        "quast/report.tsv"
    log:
        "quast/quast.log"
    params:
        output_dir = f"-o quast",
        organism = f"--{organism}" if organism != "prokaryote" else "",
        quast_params = "--labels spades_contigs,athena_scaffolds --glimmer --rna-finding" 
    threads:
        workflow.cores
    conda:
        "envs/assembly.yaml"
    shell:
        "metaquast.py --threads {threads} --pe1 {input.fastq_f} --pe2 {input.fastq_r} {params} {input.contigs} {input.scaffolds} 2> {log}"

rule BUSCO_analysis:
    input:
        "athena/athena.asm.fa"
    output:
        f"busco/short_summary.specific.{lineagedb}_odb{odb_version}.busco.txt"
    log:
        "logs/busco.log"
    params:
        method = "-m genome",
        output_folder = f"--out_path .",
        out_prefix = "-o busco",
        db_location = f"--download_path busco",
        lineage = f"-l {lineagedb}",
        metaeuk = "--metaeuk" if organism == "eukaryote" else "" 
    threads:
        workflow.cores
    conda:
        "envs/assembly.yaml"
    shell:
        """
        ( busco -f -i {input} -c {threads} {params} > {log} 2>&1 ) || touch {output}
        """

rule build_report:
    input:
        f"busco/short_summary.specific.{lineagedb}_odb{odb_version}.busco.txt",
        "quast/report.tsv"
    output:
        "reports/assembly.metrics.html"
    log:
        "logs/multiqc.log"
    params:
        options = "-n stdout --no-ai --no-version-check --force --quiet --no-data-dir",
        title = "--title \"Metassembly Metrics\""
    conda:
        "envs/qc.yaml"
    shell:
        "multiqc {params} {input} > {output} 2> {log}"

rule all:
    default_target: True
    input:
        "athena/athena.asm.fa",
        "reports/assembly.metrics.html" if not skip_reports else []
