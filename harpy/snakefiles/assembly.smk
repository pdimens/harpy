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

rule cloudspades:
    input:
        FQ_R1 = FQ1,
        FQ_R2 = FQ2
    output:
        f"{outdir}/contigs.fasta",
        f"{outdir}/scaffolds.fasta"
    params:
        outdir = outdir,
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
        "spades.py -t {threads} -m {params.mem} -k {params.k} {params.extra} --gemcode1-1 {input.FQ_R1C} --gemcode1-2 {input.FQ_R2C} -o {params.outdir} --isolate > {log}"

rule workflow_summary:
    default_target: True
    input:
        f"{outdir}/athena/athena.asm.fa"
    params:
        k_param = k_param,
        max_mem = max_mem,
        extra = extra
    run:
        with open(outdir + "/workflow/metassembly.summary", "w") as f:
            _ = f.write("The harpy assemble workflow ran using these parameters:\n\n")
            _ = f.write(f"Reads were assembled using cloudspades:\n")
            _ = f.write(f"    spades.py -t THREADS -m MEM --gemcode1-1 FQ1 --gemcode1-2 FQ2 --isolate -k {params.k_param} {params.extra}\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")