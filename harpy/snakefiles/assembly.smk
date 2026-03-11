import os

VERSION = config['Workflow']['harpy-version']
FQ1 = config["Inputs"]["fastq-r1"]
FQ2 = config["Inputs"]["fastq-r2"]
skip_reports  = config["Workflow"]["reports"]["skip"]
organism = config["Workflow"]["reports"]["organism-type"]
lineage_map = {
    "eukaryote": "eukaryota",
    "fungus": "fungi",
    "bacteria": "bacteria"
}
lineagedb = lineage_map.get(organism, "bacteria")
odb_version = 12

# SPADES
max_mem      = config["Parameters"]["spades"].get("max-memory", 'auto')
k_param      = config["Parameters"]["spades"].get("k", 10000)
spades_extra = config["Parameters"]["spades"].get("extra", "")

# ARCS
mapq       = config["Parameters"]["tigmint"].get("minimum_mapping-quality", 0)
mismatch   = config["Parameters"]["tigmint"].get("mismatch", 5)
mol_dist   = config["Parameters"]["tigmint"].get("molecule-distance", 50000)
mol_len    = config["Parameters"]["tigmint"].get("molecule-length", 2000)
span       = config["Parameters"]["tigmint"].get("span", 20)
min_align  = config["Parameters"]["arcs"].get("minimum-aligned-reads", 5)
min_contig = config["Parameters"]["arcs"].get("minimum-contig-length", 500)
seq_id     = config["Parameters"]["arcs"].get("minimum-sequence-identity", 98)
arcs_extra = config["Parameters"]["arcs"].get("extra", "")
links      = config["Parameters"]["links"].get("minimum-links", 5)

rule cloudspades:
    input:
        FQ_R1 = FQ1,
        FQ_R2 = FQ2
    output:
        "spades/contigs.fasta",
        "spades/scaffolds.fasta"
    params:
        outdir = "spades",
        k = k_param,
        mem = max_mem // 1000,
        extra = spades_extra
    log:
        "logs/assembly.log"
    conda:
        "envs/assembly.yaml"
    container:
        f"docker://pdimens/harpy:assembly_{VERSION}"
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
        temp("scaffold/interleaved.fq.gz")
    shell:
        "seqtk mergepe {input} | bgzip > {output}"

rule link_assembly:
    input:
        "spades/scaffolds.fasta",
    output:
        "scaffold/spades.fa"
    shell:  
        "ln -sr {input} {output}"

rule scaffolding:
    input:
        asm = "scaffold/spades.fa",
        reads = "scaffold/interleaved.fq.gz"
    output:
        "scaffolds.fasta"
    log:
        "logs/scaffolding.log"
    threads:
        workflow.cores
    params:
        workdir = "scaffold",
        threads = f"-j {workflow.cores}",
        draft_asm = "draft=spades",
        reads = "reads=interleaved",
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
        "envs/assembly.yaml"
    container:
        f"docker://pdimens/harpy:assembly_{VERSION}"
    shell:
        """
        arcs-make arcs-tigmint -C {params} 2> {log}
        mv {params.workdir}/spades.tigmint*.scaffolds.fa {output}
        """

rule QUAST_assessment:
    input:
        contigs = "spades/contigs.fasta",
        scaffolds = "scaffolds.fasta",
        fastq = "scaffold/interleaved.fq.gz"
    output:
        "quast/report.tsv"
    log:
        "quast/quast.log"
    params:
        output_dir = "-o quast",
        organism = f"--{organism}" if organism != "prokaryote" else "",
        quast_params = "--labels spades_contigs,arcs_scaffolds --rna-finding",
        skip_things = "--no-sv"
    threads:
        workflow.cores
    conda:
        "envs/assembly.yaml"
    container:
        f"docker://pdimens/harpy:assembly_{VERSION}"
    shell:
        "quast.py --threads {threads} --pe12 {input.fastq} {params} {input.contigs} {input.scaffolds} 2> {log}"

rule BUSCO_analysis:
    input:
        f"scaffolds.fasta"
    output:
        f"busco/short_summary.specific.{lineagedb}_odb{odb_version}.busco.txt"
    log:
        "logs/busco.log"
    params:
        #output_folder = f"--out_path {outdir}",
        out_prefix = "-o busco",
        lineage = f"-l {lineagedb}_odb{odb_version}",
        download_path = "--download_path busco",
        metaeuk = "--metaeuk" if organism == "eukaryote" else "" 
    threads:
        workflow.cores
    conda:
        "envs/assembly.yaml"
    container:
        f"docker://pdimens/harpy:assembly_{VERSION}"
    shell:
        "( busco -f -i {input} -c {threads} -m genome {params} > {log} 2>&1 ) || touch {output}"

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
        title = "--title \"Assembly Metrics\""
    conda:
        "envs/qc.yaml"
    container:
        f"docker://pdimens/harpy:qc_{VERSION}"
    shell:
        "multiqc {params} {input} > {output} 2> {log}"

rule all:
    default_target: True
    input:
        "scaffolds.fasta",
        "reports/assembly.metrics.html" if not skip_reports else []
