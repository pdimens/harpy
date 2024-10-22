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
envdir = os.path.join(os.getcwd(), ".harpy_envs")
max_mem = config["spades"]["max_memory"]
k_param = config["spades"]["k"]
metassembly = config["spades"]["assembler"]
extra = config["spades"].get("extra", "")
cloudspades = True if metassembly == "cloudspades" else False
spadesdir = f"{outdir}/{'cloudspades' if cloudspades else 'spades'}_assembly"

rule sort_by_barcode:
    input:
        fq_f = FQ1,
        fq_r = FQ2
    output:
        fq_f = temp(f"{outdir}/fastq_preproc/tmp.R1.fq"),
        fq_r = temp(f"{outdir}/fastq_preproc/tmp.R2.fq")
    params:
        barcode_tag = config["barcode_tag"].upper()
    threads:
        workflow.cores
    container:
        None
    shell:
        """
        samtools import -T "*" {input} |
        samtools sort -@ {threads} -O SAM -t {params.barcode_tag} |
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
    threads:
        workflow.cores
    resources:
        mem_mb=max_mem
    conda:
        f"{envdir}/spades.yaml"
    container:
        None
    shell:
        "metaspades.py -t {threads} -m {params.mem} -k {params.k} {params.extra} -1 {input.FQ_R1} -2 {input.FQ_R2} -o {params.outdir} --only-error-correction > {log}"

rule spades_assembly:
    input:
        fastq_R1C = outdir + "/error_correction/corrected/input.R1.fq00.0_0.cor.fastq.gz",
        fastq_R2C = outdir + "/error_correction/corrected/input.R2.fq00.0_0.cor.fastq.gz",
        fastq_UNC = outdir + "/error_correction/corrected/input.R_unpaired00.0_0.cor.fastq.gz"
    output:
        f"{outdir}/spades_assembly/contigs.fasta" 
    params:
        outdir = outdir + "/spades_assembly",
        k = k_param,
        mem = max_mem // 1000,
        extra = extra
    log:
        outdir + "/logs/spades_assembly.log"
    threads:
        workflow.cores
    resources:
        mem_mb=max_mem
    conda:
        f"{envdir}/spades.yaml"
    container:
        None
    shell:
        "metaspades.py -t {threads} -m {params.mem} -k {params.k} {params.extra} -1 {input.fastq_R1C} -2 {input.fastq_R2C} -s {input.fastq_UNC} -o {params.outdir} --only-assembler > {log}"

rule cloudspades_assembly:
    input:
        fastq_R1 = FQ1,
        fastq_R2 = FQ2
    output:
        f"{outdir}/cloudspades_assembly/contigs.fasta",
        f"{outdir}/cloudspades_assembly/scaffolds.fasta"
    params:
        outdir = spadesdir,
        k = k_param,
        mem = max_mem // 1000,
        extra = extra
    log:
        outdir + "/logs/assembly.log"
    conda:
        f"{envdir}/assembly.yaml"
    threads:
        workflow.cores
    resources:
        mem_mb=max_mem
    shell:
        "spades.py --meta -t {threads} -m {params.mem} -k {params.k} {params.extra} --gemcode1-1 {input.fastq_R1} --gemcode1-2 {input.fastq_R2} -o {params.outdir} > {log}"

rule bwa_index:
    input:
        f"{spadesdir}/contigs.fasta"
    output:
        multiext(f"{spadesdir}/contigs.fasta.", "ann", "bwt", "pac", "sa", "amb") 
    log:
        f"{outdir}/logs/bwa.index.log"
    conda:
        f"{envdir}/align.yaml"
    shell:
        "bwa index {input}"

rule bwa_align:
    input:
        multiext(f"{spadesdir}/contigs.fasta.", "ann", "bwt", "pac", "sa", "amb"),
        fastq   = collect(outdir + "/fastq_preproc/input.R{X}.fq.gz", X = [1,2]),
        contigs = f"{spadesdir}/contigs.fasta"
    output:
        f"{outdir}/reads-to-spades.bam"
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
        f"{outdir}/reads-to-spades.bam"
    output:
       f"{outdir}/reads-to-spades.bam.bai"
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
        f"{outdir}/reads-to-spades.bam.bai",
        fastq = f"{outdir}/fastq_preproc/interleaved.fq",
        bam = f"{outdir}/reads-to-spades.bam",
        contigs = f"{spadesdir}/contigs.fasta"
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
        multiext(f"{outdir}/reads-to-spades.", "bam", "bam.bai"),
        f"{outdir}/fastq_preproc/interleaved.fq",
        f"{spadesdir}/contigs.fasta",
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
        extra = extra
    run:
        summary = ["The harpy metassembly workflow ran using these parameters:"]  
        bxsort = "FASTQ inputs were sorted by their linked-read barcodes:\n"
        bxsort += "\tsamtools import -T \"*\" FQ1 FQ2 |\n"
        bxsort += f"\tsamtools sort -O SAM -t {params.bx} |\n"  
        bxsort += "\tsamtools fastq -T \"*\" -1 FQ_out1 -2 FQ_out2"  
        summary.append(bxsort)
        bxappend = "Barcoded-sorted FASTQ files had \"-1\" appended to the barcode to make them Athena-compliant:\n"  
        bxappend += f"\tsed 's/{params.bx}:Z:[^[:space:]]*/&-1/g' FASTQ | bgzip > FASTQ_OUT"  
        summary.append(bxappend)
        spades = f"Reads were assembled using {metassembly}:\n"
        if cloudspades:
            spades += f"\tspades.py -t THREADS -m {max_mem} --gemcode1-1 FQ1 --gemcode1-2 FQ2 --meta -k {k_param} {params.extra}"
        else:
            spades += f"\tmetaspades.py -t THREADS -m {max_mem} -k {k_param} {extra} -1 FQ_1 -2 FQ2 -o {spadesdir}"
        summary.append(spades)
        align = "Original input FASTQ files were aligned to the metagenome using BWA:\n"
        align += "\tbwa mem -C -p spades.contigs FQ1 FQ2 | samtools sort -O bam -"
        summary.append(align)
        interleaved = "Barcode-sorted Athena-compliant sequences were interleaved with seqtk:\n"
        interleaved += "\tseqtk mergepe FQ1 FQ2 > INTERLEAVED.FQ"
        summary.append(interleaved)
        athena = "Athena ran with the config file Harpy built from the files created from the previous steps:\n"
        athena += "\tathena-meta --config athena.config"
        summary.append(athena)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['workflow_call']}"
        summary.append(sm)
        with open(outdir + "/workflow/metassembly.summary", "w") as f:  
            f.write("\n\n".join(summary))
