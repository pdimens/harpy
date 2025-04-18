containerized: "docker://pdimens/harpy:latest"

import os
import gzip
import shutil
import logging
from pathlib import Path
from itertools import product

envdir   = os.path.join(os.getcwd(), "workflow", "envs")
gen_hap1 = config["inputs"]["genome_hap1"]
gen_hap2 = config["inputs"]["genome_hap2"]
barcode_file = config["barcodes"]["file"]
barcode_len = config["barcodes"]["length"]
merge_haplotypes = config["merge_haplotypes"]
genodict = {"0": gen_hap1, "1": gen_hap2}

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake_log"])
    logger.addHandler(logfile_handler)
onsuccess:
    shutil.rmtree('_Inline', ignore_errors=True)
onerror:
    shutil.rmtree('_Inline', ignore_errors=True)
wildcard_constraints:
    hap = r"[01]"

rule barcode_keymap:
    input:
        bc = barcode_file
    output:
        bc = "barcodes.key.gz"
    run:
        bc_range = [f"{i}".zfill(2) for i in range(1,97)]
        bc_generator = product("A", bc_range, "C", bc_range, "B", bc_range, "D", bc_range)
        with open(input.bc, "r") as bc_in, gzip.open(output.bc, "wb") as bc_out:
            for nuc_barcode in bc_in:
                haptag = "".join(next(bc_generator))
                bc_out.write((nuc_barcode.rstrip() + "\t" + haptag + "\n").encode("utf-8"))

rule link_genome:
    input:
        hap = lambda wc: genodict[wc.get("hap")]
    output: 
        fa = "workflow/input/hap.{hap}.fasta"
    run:
        try:
            with gzip.open(input.hap, 'rb') as gzip_file, open(output.fa, 'wb') as output_file:
                shutil.copyfileobj(gzip_file, output_file)
        except gzip.BadGzipFile:
            if not (Path(output.fa).is_symlink() or Path(output.fa).exists()):
                Path(output.fa).symlink_to(Path(input.hap).resolve()) 

rule index_genome:
    input:
        "workflow/input/hap.{hap}.fasta"
    output: 
        "workflow/input/hap.{hap}.fasta.fai"
    container:
        None
    shell:
        "samtools faidx --fai-idx {output} {input}"

if barcode_file == "workflow/input/haplotag_barcodes.txt":
    rule create_barcodes:
        output:
            "workflow/input/haplotag_barcodes.txt"
        container:
            None
        shell:
            "haplotag_barcodes.py > {output}"

rule simulate_reads:
    input:
        "workflow/input/hap.{hap}.fasta"
    output:
        temp(multiext("dwgsim/sim_reads.{hap}.12", ".bwa.read1.fastq.gz" ,".bwa.read2.fastq.gz", ".mutations.txt", ".mutations.vcf"))
    log:
        "logs/dwgsim.hap.{hap}.log"
    params:
        readpairs = int(config["read_pairs"] * 500000),
        outerdist = f"-d {config['outer_distance']}",
        static = "-e 0.0001,0.0016 -E 0.0001,0.0016 -1 135 -2 151 -H -y 0 -S 0 -c 0 -R 0 -o 1 -F 0 -m /dev/null",
        distsd = f"-s {config['distance_sd']}",
        mutationrate = f"-r {config['mutation_rate']}",
        input_file = lambda wc: "workflow/input/hap." + wc.get("hap") + ".fasta",
        prefix = lambda wc: "dwgsim/sim_reads." + wc.get("hap") + ".12"
    conda:
        f"{envdir}/simulations.yaml"
    shell:
        "dwgsim -N {params} 2> {log}"

rule interleave_reads:
    input:
        collect("dwgsim/sim_reads.{{hap}}.12.bwa.read{rd}.fastq.gz", rd = [1,2]) 
    output:
        "dwgsim/sim_reads.{hap}.12.fastq"
    container:
        None
    shell:
        "seqtk mergepe {input} > {output}"

rule create_molecules:
    input:
        hap_reads = collect("dwgsim/sim_reads.{hap}.12.fastq"   , hap = [0,1]),
        fasta_fai = collect("workflow/input/hap.{hap}.fasta.fai", hap = [0,1]),
        barcodes = barcode_file
    output:
        temp(collect("linked_molecules/lrsim.{hap}.fp"      , hap = [0,1])),
        temp(collect("linked_molecules/lrsim.{hap}.manifest", hap = [0,1]))
    log:
        "logs/linked_molecules.log"
    params:
        haplosim = "workflow/scripts/HaploSim.pl",
        reads_in = f"-a dwgsim/sim_reads.0.12.fastq,dwgsim/sim_reads.1.12.fastq",
        fai_in   = f"-g workflow/input/hap.0.fasta.fai,workflow/input/hap.1.fasta.fai",
        bccodes  = f"-b {barcode_file}",
        proj_dir = f"-p linked_molecules/lrsim",
        outdist  = f"-i {config['outer_distance']}",
        dist_sd  = f"-s {config['distance_sd']}",
        n_pairs  = f"-x {config['read_pairs']}",
        mol_len  = f"-f {config['molecule_length']}",
        parts    = f"-t {config['partitions']}",
        mols_per = f"-m {config['molecules_per_partition']}",
        bc_len   = f"-l {barcode_len}"
    threads:
        workflow.cores
    conda:
        f"{envdir}/simulations.yaml"
    shell: 
        "perl {params} -z {threads} -d 2 2> {log}"

rule sort_molecules:
    input:
        "linked_molecules/lrsim.{hap}.manifest"
    output:
        "linked_molecules/lrsim.{hap}.sort.manifest"
    conda:
        f"{envdir}/simulations.yaml"
    shell:
        "msort -kn1 {input} > {output}"

rule create_linked_reads:
    input:
        manifest = "linked_molecules/lrsim.{hap}.sort.manifest",
        dwg_hap = "dwgsim/sim_reads.{hap}.12.fastq"
    output:
        collect("multiplex/sim_hap{{hap}}_multiplex_R{FR}_001.fastq.gz", FR = [1,2])
    log:
        "logs/create_linkedreads.hap{hap}.log"
    params:
        lambda wc: f"multiplex/sim_hap{wc.get('hap')}_multiplex"
    container:
        None
    shell:
        "extractReads {input} {params} 2> {log}"

rule demultiplex_barcodes:
    input:
        fw = "multiplex/sim_hap{hap}_multiplex_R1_001.fastq.gz",
        rv = "multiplex/sim_hap{hap}_multiplex_R2_001.fastq.gz",
        barcodes = "barcodes.key.gz"
    output:
        fw = "sim_hap{hap}.R1.fq.gz",
        rv = "sim_hap{hap}.R2.fq.gz"
    log:
        "logs/sim_hap{hap}.demultiplex"
    params:
        lambda wc: f"sim_hap{wc.get('hap')}"
    container:
        None
    shell:
        "inline_to_haplotag.py -b {input.barcodes} -p {params} {input.fw} {input.rv} 2> {log}"

rule combine_haplotypes:
    input:
        hap1 = "sim_hap0.R{FR}.fq.gz",
        hap2 = "sim_hap1.R{FR}.fq.gz",
    output:
        "sim_hap.R{FR}.fq.gz"
    container:
        None
    shell:
        "cat {input} > {output}"

rule workflow_summary:
    default_target: True
    input:
        collect("sim_hap{hap}.R{fw}.fq.gz", hap = [0,1], fw = [1,2]),
        collect("sim_hap.R{FR}.fq.gz", FR = [1,2]) if merge_haplotypes else []
    params:
        lrsproj_dir = os.getcwd(),
        lrsoutdist  = config["outer_distance"],
        lrsdist_sd  = config["distance_sd"],
        lrsn_pairs  = config["read_pairs"],
        lrsmol_len  = config["molecule_length"],
        lrsparts    = config["partitions"],
        lrsmols_per = config["molecules_per_partition"],
        lrbc_len   = barcode_len,
        lrsstatic = "-o 1 -d 2",
        dwgreadpairs = int(config["read_pairs"] * 500000),
        dwgouterdist = config["outer_distance"],
        dwgdistsd = config["distance_sd"],
        dwgmutationrate = config["mutation_rate"],
        dwgprefix = "dwgsim/sim_reads.hap.12"
    run:
        summary = ["The harpy simulate linkedreas workflow ran using these parameters:"]
        summary.append(f"Genome haplotype 1: {gen_hap1}")
        summary.append(f"Genome haplotype 2: {gen_hap2}")
        summary.append(f"Barcode file: {barcode_file}")
        dwgsim = "Reads were simulated from the provided genomes using:\n"
        dwgsim += f"\tdwgsim -N {params.dwgreadpairs} -e 0.0001,0.0016 -E 0.0001,0.0016 -d {params.dwgouterdist} -s {params.dwgdistsd} -1 135 -2 151 -H -y 0 -S 0 -c 0 -R 0 -r {params.dwgmutationrate} -F 0 -o 1 -m /dev/null GENO PREFIX"
        summary.append(dwgsim)
        haplosim = "HaploSim (Harpy's fork of LRSIM) was used with these parameters:\n"
        haplosim += f"\tHaploSim.pl -g genome1,genome2 -a dwgsimreads1,dwgsimreads2 -l {params.lrbc_len} -p {params.lrsproj_dir}/linked_molecules/lrsim -b BARCODES -i {params.lrsoutdist} -s {params.lrsdist_sd} -x {params.lrsn_pairs} -f {params.lrsmol_len} -t {params.lrsparts} -m {params.lrsmols_per} -z THREADS {params.lrsstatic}"
        summary.append(haplosim)
        bxconvert = "Inline barcodes were converted in haplotag BX:Z tags using:\n"
        bxconvert += "\tinline_to_haplotag.py -b <barcodes.txt> -p <prefix> forward.fq.gz reverse.fq.gz"
        summary.append(bxconvert)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['snakemake_command']}"
        summary.append(sm)
        with open("workflow/simulate.reads.summary", "w") as f:
            f.write("\n\n".join(summary))
