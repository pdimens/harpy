containerized: "docker://pdimens/harpy:latest"

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
# SPADES
max_mem      = config["spades"]["max_memory"]
k_param      = config["spades"]["k"]
spades_extra = config["spades"].get("extra", "")
# ARCS
mapq       = config["tigmint"]["minimum_mapping_quality"]
mismatch   = config["tigmint"]["mismatch"]
mol_dist   = config["tigmint"]["molecule_distance"]
mol_len    = config["tigmint"]["molecule_length"]
span       = config["tigmint"]["span"]
min_align  = config["arcs"]["minimum_aligned_reads"]
min_contig = config["arcs"]["minimum_contig_length"]
seq_id     = config["arcs"]["minimum_sequence_identity"]
arcs_extra = config["arcs"].get("extra", "")
links      = config["links"]["minimum_links"]

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
        extra = spades_extra
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
    container:
        None
    shell:
        "seqtk mergepe {input} > {output}"

rule link_assembly:
    input:
        f"{outdir}/spades/scaffolds.fasta",
    output:
        f"{outdir}/scaffold/spades.fa"
    container:
        None
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
        min_mapq = f"mapq={mapq}",
        max_mismatch = f"nm={mismatch}",
        moldist = f"dist={mol_dist}",
        min_length = f"minsize={mol_len}",
        span = f"span={span}",
        min_perbarcod = f"c={min_align}",
        min_contig = f"z={min_contig}",
        min_seqid = f"s={seq_id}",
        min_links = f"l={links}",
        extra = arcs_extra
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
        extra = extra,
        workdir = f"-C {outdir}/scaffold",
        threads = f"-j {workflow.cores}",
        draft_asm = f"draft=spades",
        reads = f"reads=interleaved",
        bwa_threads = f"t={workflow.cores}",
        min_mapq = f"mapq={mapq}",
        max_mismatch = f"nm={mismatch}",
        moldist = f"dist={mol_dist}",
        min_length = f"minsize={mol_len}",
        span = f"span={span}",
        min_perbarcod = f"c={min_align}",
        min_contig = f"z={min_contig}",
        min_seqid = f"s={seq_id}",
        min_links = f"l={links}",
        extra = arcs_extra
    run:
        summary_template = f"""
The harpy assemble workflow ran using these parameters:

Reads were assembled using cloudspades:
    spades.py -t THREADS -m {params.max_mem} --gemcode1-1 FQ1 --gemcode1-2 FQ2 --isolate -k {params.k_param} {params.extra}

The draft assembly was error corrected and scaffolded with ARCS:
    arcs-make arcs {params[3:]}

The Snakemake workflow was called via command line:
    {config["workflow_call"]}
"""
        with open(outdir + "/workflow/metassembly.summary", "w") as f:
            f.write(summary_template)