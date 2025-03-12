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
envdir = os.path.join(os.getcwd(), outdir, "workflow", "envs")
skip_reports  = config["reports"]["skip"]
organism = config["reports"]["organism_type"]
lineage_map = {
    "eukaryote": "eukaryota",
    "fungus": "fungi",
    "bacteria": "bacteria"
}
lineagedb = lineage_map.get(organism, "bacteria")

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
        temp(f"{outdir}/scaffold/interleaved.fq.gz")
    container:
        None
    shell:
        "seqtk mergepe {input} | bgzip > {output}"

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
    input:
        asm = f"{outdir}/scaffold/spades.fa",
        reads = f"{outdir}/scaffold/interleaved.fq.gz"
    output:
        f"{outdir}/scaffolds.fasta"
    log:
        outdir + "/logs/scaffolding.log"
    threads:
        workflow.cores
    params:
        workdir = f"{outdir}/scaffold",
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
        prefix = "base_name=scaffolds",
        extra = arcs_extra
    conda:
        f"{envdir}/assembly.yaml"
    shell:
        """
        arcs-make arcs-tigmint -C {params} 2> {log}
        mv {params.workdir}/spades.tigmint*.scaffolds.fa {output}
        """

rule QUAST_assessment:
    input:
        contigs = f"{outdir}/spades/contigs.fasta",
        scaffolds = f"{outdir}/scaffolds.fasta",
        fastq = f"{outdir}/scaffold/interleaved.fq.gz"
    output:
        f"{outdir}/quast/report.tsv"
    log:
        f"{outdir}/quast/quast.log"
    params:
        output_dir = f"-o {outdir}/quast",
        organism = f"--{organism}" if organism != "prokaryote" else "",
        quast_params = "--labels spades_contigs,arcs_scaffolds --rna-finding",
        skip_things = "--no-sv"
    threads:
        workflow.cores
    conda:
        f"{envdir}/assembly.yaml"
    shell:
        "quast.py --threads {threads} --pe12 {input.fastq} {params} {input.contigs} {input.scaffolds} 2> {log}"

rule BUSCO_analysis:
    input:
        f"{outdir}/scaffolds.fasta"
    output:
        f"{outdir}/busco/short_summary.specific.{lineagedb}_odb12.busco.txt"
    log:
        f"{outdir}/logs/busco.log"
    params:
        output_folder = outdir,
        out_prefix = "-o busco",
        lineage = f"-l {lineagedb}_odb12",
        download_path = f"--download_path {outdir}/busco",
        metaeuk = "--metaeuk" if organism == "eukaryote" else "" 
    threads:
        workflow.cores
    conda:
        f"{envdir}/assembly.yaml"
    shell:
        "( busco -f -i {input} -c {threads} -m genome --out_path {params} > {log} 2>&1 ) || touch {output}"

rule build_report:
    input:
        f"{outdir}/busco/short_summary.specific.{lineagedb}_odb12.busco.txt",
        f"{outdir}/quast/report.tsv"
    output:
        f"{outdir}/reports/assembly.metrics.html"
    params:
        options = "--no-version-check --force --quiet --no-data-dir",
        title = "--title \"Assembly Metrics\""
    conda:
        f"{envdir}/qc.yaml"
    shell:
        "multiqc {input} {params} --filename {output}"

rule workflow_summary:
    default_target: True
    input:
        f"{outdir}/scaffolds.fasta",
        f"{outdir}/reports/assembly.metrics.html" if not skip_reports else [],
    params:
        k_param = k_param,
        max_mem = max_mem // 1000,
        spades_extra = spades_extra,
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
        arcs_extra = arcs_extra
    run:
        summary = ["The harpy assemble workflow ran using these parameters:"]
        spades = "Reads were assembled using cloudspades:\n"
        spades += f"\tspades.py -t THREADS -m {params.max_mem} --gemcode1-1 FQ1 --gemcode1-2 FQ2 --isolate -k {params.k_param} {params.spades_extra}"
        summary.append(spades)
        arcs = "The draft assembly was error corrected and scaffolded with Tigmint/ARCS/LINKS:\n"
        arcs += f"\tarcs-make arcs-tigmint {" ".join(params[3:])}"
        summary.append(arcs)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['workflow_call']}"
        summary.append(sm)
        with open(outdir + "/workflow/assembly.summary", "w") as f:
            f.write("\n\n".join(summary))