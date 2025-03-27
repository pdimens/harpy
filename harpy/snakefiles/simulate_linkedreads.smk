containerized: "docker://pdimens/harpy:latest"

import os
import gzip
import shutil
import logging
from pathlib import Path
from itertools import product

outdir   = config["output_directory"]
envdir   = os.path.join(os.getcwd(), outdir, "workflow", "envs")
gen_hap1 = config["inputs"]["genome_hap1"]
gen_hap2 = config["inputs"]["genome_hap2"]
barcode_file = config["barcodes"]["file"]
barcode_len = config["barcodes"]["length"]
merge_haplotypes = config["merge_haplotypes"]
genodict = {"0": gen_hap1, "1": gen_hap2}

onstart:
    logger.addHandler(logging.FileHandler(config["snakemake_log"]))
onsuccess:
    os.remove(logger.logfile)
    shutil.rmtree('_Inline', ignore_errors=True)
onerror:
    os.remove(logger.logfile)
    shutil.rmtree('_Inline', ignore_errors=True)
wildcard_constraints:
    hap = r"[01]"

rule barcode_keymap:
    input:
        barcode_file
    output:
        outdir + "/barcodes.key.gz"
    run:
        bc_range = [f"{i}".zfill(2) for i in range(1,97)]
        bc_generator = product("A", bc_range, "C", bc_range, "B", bc_range, "D", bc_range)
        with open(input[0], "r") as bc_in, gzip.open(output[0], "wb") as bc_out:
            for nuc_barcode in bc_in:
                haptag = "".join(next(bc_generator))
                bc_out.write((nuc_barcode.rstrip() + "\t" + haptag + "\n").encode("utf-8"))

rule link_genome:
    input:
        lambda wc: genodict[wc.get("hap")]
    output: 
        outdir + "/workflow/input/hap.{hap}.fasta"
    run:
        if input[0].lower().endswith("gz"):
            with gzip.open(input[0], 'rb') as gzip_file, open(output[0], 'wb') as output_file:
                shutil.copyfileobj(gzip_file, output_file)
        else:
            if not (Path(output[0]).is_symlink() or Path(output[0]).exists()):
                Path(output[0]).symlink_to(Path(input[0]).resolve()) 

rule index_genome:
    input:
        outdir + "/workflow/input/hap.{hap}.fasta"
    output: 
        outdir + "/workflow/input/hap.{hap}.fasta.fai"
    container:
        None
    shell:
        "samtools faidx --fai-idx {output} {input}"

if barcode_file == f"{outdir}/workflow/input/haplotag_barcodes.txt":
    rule create_barcodes:
        output:
            f"{outdir}/workflow/input/haplotag_barcodes.txt"
        container:
            None
        shell:
            "haplotag_barcodes.py > {output}"

rule simulate_reads:
    input:
        outdir + "/workflow/input/hap.{hap}.fasta"
    output:
        temp(multiext(outdir + "/dwgsim/sim_reads.{hap}.12", ".bwa.read1.fastq.gz" ,".bwa.read2.fastq.gz", ".mutations.txt", ".mutations.vcf"))
    log:
        outdir + "/logs/dwgsim.hap.{hap}.log"
    params:
        readpairs = int(config["read_pairs"] * 500000),
        outerdist = config["outer_distance"],
        distsd = config["distance_sd"],
        mutationrate = config["mutation_rate"],
        prefix = lambda wc: outdir + "/dwgsim/sim_reads." + wc.get("hap") + ".12"
    conda:
        f"{envdir}/simulations.yaml"
    shell:
        """
        dwgsim -N {params.readpairs} -e 0.0001,0.0016 -E 0.0001,0.0016 -d {params.outerdist} -s {params.distsd} -1 135 -2 151 -H -y 0 -S 0 -c 0 -R 0 -o 1 -r {params.mutationrate} -F 0 -m /dev/null {input} {params.prefix} 2> {log}
        """

rule interleave_reads:
    input:
        collect(outdir + "/dwgsim/sim_reads.{{hap}}.12.bwa.read{rd}.fastq.gz", rd = [1,2]) 
    output:
        outdir + "/dwgsim/sim_reads.{hap}.12.fastq"
    container:
        None
    shell:
        "seqtk mergepe {input} > {output}"

rule create_molecules:
    input:
        hap_reads = collect(outdir + "/dwgsim/sim_reads.{hap}.12.fastq"   , hap = [0,1]),
        fasta_fai = collect(outdir + "/workflow/input/hap.{hap}.fasta.fai", hap = [0,1]),
        barcodes = barcode_file
    output:
        temp(collect(outdir + "/linked_molecules/lrsim.{hap}.fp"      , hap = [0,1])),
        temp(collect(outdir + "/linked_molecules/lrsim.{hap}.manifest", hap = [0,1]))
    log:
        f"{outdir}/logs/linked_molecules.log"
    params:
        haplosim = f"{outdir}/workflow/scripts/HaploSim.pl",
        reads_in = f"-a {outdir}/dwgsim/sim_reads.0.12.fastq,{outdir}/dwgsim/sim_reads.1.12.fastq",
        fai_in   = f"-g {outdir}/workflow/input/hap.0.fasta.fai,{outdir}/workflow/input/hap.1.fasta.fai",
        bccodes  = f"-b {barcode_file}",
        proj_dir = f"-p {outdir}/linked_molecules/lrsim",
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
        outdir + "/linked_molecules/lrsim.{hap}.manifest"
    output:
        outdir + "/linked_molecules/lrsim.{hap}.sort.manifest"
    conda:
        f"{envdir}/simulations.yaml"
    shell:
        "msort -kn1 {input} > {output}"

rule create_linked_reads:
    input:
        manifest = outdir + "/linked_molecules/lrsim.{hap}.sort.manifest",
        dwg_hap = outdir + "/dwgsim/sim_reads.{hap}.12.fastq"
    output:
        collect(outdir + "/multiplex/sim_hap{{hap}}_multiplex_R{FR}_001.fastq.gz", FR = [1,2])
    log:
        outdir + "/logs/create_linkedreads.hap{hap}.log"
    params:
        lambda wc: f"{outdir}/multiplex/sim_hap{wc.get('hap')}_multiplex"
    container:
        None
    shell:
        "extractReads {input} {params} 2> {log}"

rule demultiplex_barcodes:
    input:
        fw = outdir + "/multiplex/sim_hap{hap}_multiplex_R1_001.fastq.gz",
        rv = outdir + "/multiplex/sim_hap{hap}_multiplex_R2_001.fastq.gz",
        barcodes = outdir + "/barcodes.key.gz"
    output:
        fw = outdir + "/sim_hap{hap}.R1.fq.gz",
        rv = outdir + "/sim_hap{hap}.R2.fq.gz"
    log:
        outdir + "/logs/sim_hap{hap}.demultiplex"
    params:
        lambda wc: f"{outdir}/sim_hap{wc.get('hap')}"
    container:
        None
    shell:
        "inline_to_haplotag.py -b {input.barcodes} -p {params} {input.fw} {input.rv} 2> {log}"

rule combine_haplotypes:
    input:
        hap1 = outdir + "/sim_hap0.R{FR}.fq.gz",
        hap2 = outdir + "/sim_hap1.R{FR}.fq.gz",
    output:
        outdir + "/sim_hap.R{FR}.fq.gz"
    container:
        None
    shell:
        "cat {input} > {output}"

rule workflow_summary:
    default_target: True
    input:
        collect(outdir + "/sim_hap{hap}.R{fw}.fq.gz", hap = [0,1], fw = [1,2]),
        collect("sim_hap.R{FR}.fq.gz", FR = [1,2]) if merge_haplotypes else []
    params:
        lrsproj_dir = f"{outdir}",
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
        dwgprefix = outdir + "/dwgsim/sim_reads.hap.12"
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
        sm += f"\t{config['workflow_call']}"
        summary.append(sm)
        with open(outdir + "/workflow/simulate.reads.summary", "w") as f:
            f.write("\n\n".join(summary))
