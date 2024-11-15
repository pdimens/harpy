containerized: "docker://pdimens/harpy:latest"

import os
import gzip
import shutil
import logging
from pathlib import Path

onstart:
    logger.logger.addHandler(logging.FileHandler(config["snakemake_log"]))
onsuccess:
    os.remove(logger.logfile)
onerror:
    os.remove(logger.logfile)

outdir   = config["output_directory"]
gen_hap1 = config["inputs"]["genome_hap1"]
gen_hap2 = config["inputs"]["genome_hap2"]
barcode_file = config["barcodes"]["file"]
barcode_len = config["barcodes"]["length"]
envdir   = os.path.join(os.getcwd(), ".harpy_envs")
genodict = {"0": gen_hap1, "1": gen_hap2}

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
        temp(multiext(outdir + "/dwgsim_simulated/dwgsim.{hap}.12", ".bwa.read1.fastq.gz" ,".bwa.read2.fastq.gz", ".mutations.txt", ".mutations.vcf"))
    log:
        outdir + "/logs/dwgsim.hap.{hap}.log"
    params:
        readpairs = int(config["read_pairs"] * 500000),
        outerdist = config["outer_distance"],
        distsd = config["distance_sd"],
        mutationrate = config["mutation_rate"],
        prefix = lambda wc: outdir + "/dwgsim_simulated/dwgsim." + wc.get("hap") + ".12"
    conda:
        f"{envdir}/simulations.yaml"
    shell:
        """
        dwgsim -N {params.readpairs} -e 0.0001,0.0016 -E 0.0001,0.0016 -d {params.outerdist} -s {params.distsd} -1 135 -2 151 -H -y 0 -S 0 -c 0 -R 0 -o 1 -r {params.mutationrate} -F 0 -m /dev/null {input} {params.prefix} 2> {log}
        """

rule interleave_reads:
    input:
        collect(outdir + "/dwgsim_simulated/dwgsim.{{hap}}.12.bwa.read{rd}.fastq.gz", rd = [1,2]) 
    output:
        outdir + "/dwgsim_simulated/dwgsim.{hap}.12.fastq"
    container:
        None
    shell:
        "seqtk mergepe {input} > {output}"

rule make_linked_reads:
    input:
        hap1 = f"{outdir}/dwgsim_simulated/dwgsim.0.12.fastq",
        hap2 = f"{outdir}/dwgsim_simulated/dwgsim.1.12.fastq",
        fai1 = outdir + "/workflow/input/hap.0.fasta.fai",
        fai2 = outdir + "/workflow/input/hap.1.fasta.fai",
        barcodes = barcode_file
    output:
        collect(outdir + "/lrsim/sim.{hap}.{ext}", hap = [0,1], ext = ["fp", "manifest"]),
        temp(f"{outdir}/lrsim/.status")
    log:
        f"{outdir}/logs/LRSIM.log"
    params:
        lrsim = f"{outdir}/workflow/scripts/LRSIM_harpy.pl",
        infiles = f"-g {outdir}/dwgsim_simulated/dwgsim.0.12.fastq,{outdir}/dwgsim_simulated/dwgsim.1.12.fastq",
        inbarcodes = f"-b {barcode_file}",
        proj_dir = f"-p {outdir}/lrsim/sim",
        prefix = f"-r {outdir}",
        outdist  = f"-i {config['outer_distance']}",
        dist_sd  = f"-s {config['distance_sd']}",
        n_pairs  = f"-x {config['read_pairs']}",
        mol_len  = f"-f {config['molecule_length']}",
        parts    = f"-t {config['partitions']}",
        mols_per = f"-m {config['molecules_per_partition']}",
        bc_len   = f"-l {barcode_len}",
        static   = "-o 1 -d 2 -u 4"
    threads:
        workflow.cores
    conda:
        f"{envdir}/simulations.yaml"
    shell: 
        "perl {params} -z {threads} 2> {log}"

rule sort_manifest:
    input:
        outdir + "/lrsim/sim.{hap}.manifest"
    output:
        outdir + "/lrsim/sim.{hap}.sort.manifest"
    conda:
        f"{envdir}/simulations.yaml"
    shell:
        "msort -kn1 {input} > {output}"

rule extract_linked_reads:
    input:
        manifest = outdir + "/lrsim/sim.{hap}.sort.manifest",
        dwg_hap = outdir + "/dwgsim_simulated/dwgsim.{hap}.12.fastq"
    output:
        outdir + "/10X/sim_hap{hap}_10x_R1_001.fastq.gz",
        outdir + "/10X/sim_hap{hap}_10x_R2_001.fastq.gz"
    log:
        outdir + "/logs/extract_linkedreads.hap{hap}.log"
    params:
        lambda wc: f"""{outdir}/10X/sim_hap{wc.get("hap")}_10x"""
    container:
        None
    shell:
        "extractReads {input} {params} 2> {log}"

rule demultiplex_barcodes:
    input:
        fw = outdir + "/10X/sim_hap{hap}_10x_R1_001.fastq.gz",
        rv = outdir + "/10X/sim_hap{hap}_10x_R2_001.fastq.gz",
        barcodes = barcode_file
    output:
        fw = outdir + "/sim_hap{hap}_haplotag.R1.fq.gz",
        rv = outdir + "/sim_hap{hap}_haplotag.R2.fq.gz"
    log:
        outdir + "/sim_hap/{hap}_haplotag.barcodes"
    params:
        outdir = outdir,
        bc_len = barcode_len
    container:
        None
    shell:
        "inline_to_haplotag.py -l {params.bc_len} -f {input.fw} -r {input.rv} -b {input.barcodes} -p {params.outdir}/sim_hap{wildcards.hap}_haplotag"

rule workflow_summary:
    default_target: True
    input:
        collect(outdir + "/sim_hap{hap}_haplotag.R{fw}.fq.gz", hap = [0,1], fw = [1,2])
    params:
        lrsproj_dir = f"{outdir}",
        lrsoutdist  = config["outer_distance"],
        lrsdist_sd  = config["distance_sd"],
        lrsn_pairs  = config["read_pairs"],
        lrsmol_len  = config["molecule_length"],
        lrsparts    = config["partitions"],
        lrsmols_per = config["molecules_per_partition"],
        lrbc_len   = barcode_len,
        lrsstatic = "-o 1 -d 2 -u 4",
        dwgreadpairs = int(config["read_pairs"] * 500000),
        dwgouterdist = config["outer_distance"],
        dwgdistsd = config["distance_sd"],
        dwgmutationrate = config["mutation_rate"],
        dwgprefix = outdir + "/dwgsim_simulated/dwgsim.hap.12"
    run:
        summary = ["The harpy simulate linkedreas workflow ran using these parameters:"]
        summary.append(f"Genome haplotype 1: {gen_hap1}")
        summary.append(f"Genome haplotype 2: {gen_hap2}")
        summary.append(f"Barcode file: {barcode_file}")
        dwgsim = "Reads were simulated from the provided genomes using:\n"
        dwgsim += f"\tdwgsim -N {params.dwgreadpairs} -e 0.0001,0.0016 -E 0.0001,0.0016 -d {params.dwgouterdist} -s {params.dwgdistsd} -1 135 -2 151 -H -y 0 -S 0 -c 0 -R 0 -r {params.dwgmutationrate} -F 0 -o 1 -m /dev/null GENO PREFIX"
        summary.append(dwgsim)
        lrsim = "LRSIM was started from step 3 (-u 3) with these parameters:\n"
        lrsim += f"\tLRSIM_harpy.pl -g genome1,genome2 -l {params.lrbc_len} -p {params.lrsproj_dir}/lrsim/sim -b BARCODES -r {params.lrsproj_dir} -i {params.lrsoutdist} -s {params.lrsdist_sd} -x {params.lrsn_pairs} -f {params.lrsmol_len} -t {params.lrsparts} -m {params.lrsmols_per} -z THREADS {params.lrsstatic}"
        summary.append(lrsim)
        bxconvert = "Inline barcodes were converted in haplotag BX:Z tags using:\n"
        bxconvert += "\tinline_to_haplotag.py"
        summary.append(bxconvert)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['workflow_call']}"
        summary.append(sm)
        with open(outdir + "/workflow/simulate.reads.summary", "w") as f:
            f.write("\n\n".join(summary))
