import os
import logging

onstart:
    logger.logger.addHandler(logging.FileHandler(config["snakemake_log"]))
onsuccess:
    os.remove(logger.logfile)
onerror:
    os.remove(logger.logfile)

FQ1 = config["inputs"]["fastq_r1"]
FQ2 = config["inputs"]["fastq_r2"]
outdir = config["output_directory"]
envdir = os.getcwd() + "/.harpy_envs"
max_mem = config["spades"]["max_memory"]
k_param = config["spades"]["k"]
extra = config["spades"].get("extra", "") 

rule cloudspades:
    input:
        FQ_R1 = FQ1,
        FQ_R2 = FQ2
    output:
        f"{outdir}/spades/contigs.fasta",
        f"{outdir}/spades/scaffolds.fasta"
    params:
        outdir = f"{outdir}/spades",
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
        "spades.py -t {threads} -m {params.mem} -k {params.k} {params.extra} --gemcode1-1 {input.FQ_R1} --gemcode1-2 {input.FQ_R2} -o {params.outdir} --isolate > {log}"

rule interleave_fastq:
    input:
        FQ1,
        FQ2
    output:
        temp(f"{outdir}/scaffold/interleaved.fq")
    shell:
        "seqtk mergepe {input} > {output}"

rule link_assembly:
    input:
        f"{outdir}/spades/scaffolds.fasta",
    output:
        f"{outdir}/scaffold/spades.fa"
    shell:  
        "ln -sr {input} {output}"

rule scaffolding:
    default_target: True
    input:
        asm = f"{outdir}/scaffold/spades.fa",
        reads = f"{outdir}/scaffold/interleaved.fq"
#    output:
#        f"{outdir}/scaffold/scaffolds.fasta.tigmint.fa"
    threads:
        workflow.cores
    params:
        workdir = f"-C {outdir}/scaffold",
        threads = f"-j {workflow.cores}",
        draft_asm = f"draft=spades",
        reads = f"reads=interleaved",
        bwa_threads = f"t={workflow.cores}",
        span = f"span={span}",
        moldist = f"dist={moldist}",
        max_mismatch = f"nm={max_mismatch}",
        min_length = f"minsize={min_length}",
        min_mapq = f"mapq={min_mapq}",
        min_contig = f"z={min_contig}",
        min_perbarcod = f"c={min_aligned_pairs}",
        min_seqid = f"s={min_seq_id}",
        min_links = f"l={min_links}"
    conda:
        f"{envdir}/assembly.yaml"
    shell:
        "arcs-make arcs {params}"

rule workflow_summary:
    #default_target: True
    input:
        f"{outdir}/athena/athena.asm.fa"
    params:
        k_param = k_param,
        max_mem = max_mem // 1000,
        extra = extra
    run:
        summary_template = f"""
The harpy assemble workflow ran using these parameters:

Reads were assembled using cloudspades:
    spades.py -t THREADS -m {params.max_mem} --gemcode1-1 FQ1 --gemcode1-2 FQ2 --isolate -k {params.k_param} {params.extra}

The Snakemake workflow was called via command line:
    {config["workflow_call"]}
"""
        with open(outdir + "/workflow/metassembly.summary", "w") as f:
            f.write(summary_template)