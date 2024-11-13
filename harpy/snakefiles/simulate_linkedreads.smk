containerized: "docker://pdimens/harpy:latest"

import os
import gzip
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
barcodes = config["inputs"].get("barcodes", None)
envdir   = os.path.join(os.getcwd(), ".harpy_envs")
barcodefile = barcodes if barcodes else f"{outdir}/workflow/input/haplotag_barcodes.txt"

rule link_1st_geno:
    input:
        gen_hap1
    output: 
        f"{outdir}/workflow/input/hap.0.fasta"
    run:
        if input[0].lower().endswith("gz"):
            with open(input[0], 'rb') as inf, open(output[0], 'w', encoding='utf8') as outf:
                decom_str = gzip.decompress(inf.read()).decode('utf-8')
                outf.write(decom_str)
        else:
            if not (Path(output[0]).is_symlink() or Path(output[0]).exists()):
                Path(output[0]).symlink_to(Path(input[0]).resolve()) 

use rule link_1st_geno as link_2nd_geno with:
    input:
        gen_hap2
    output: 
        f"{outdir}/workflow/input/hap.1.fasta"

rule index_genome:
    input:
        outdir + "/workflow/input/hap.{hap}.fasta"
    output: 
        outdir + "/workflow/input/hap.{hap}.fasta.fai"
    container:
        None
    shell:
        "samtools faidx --fai-idx {output} {input}"

if not barcodes:
    rule create_barcodes:
        output:
            barcodefile
        container:
            None
        shell:
            "haplotag_barcodes.py > {output}"

rule create_reads:
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

rule interleave_dwgsim:
    input:
        collect(outdir + "/dwgsim_simulated/dwgsim.{{hap}}.12.bwa.read{rd}.fastq.gz", rd = [1,2]) 
    output:
        outdir + "/dwgsim_simulated/dwgsim.{hap}.12.fastq"
    container:
        None
    shell:
        "seqtk mergepe {input} > {output}"

rule lrsim:
    input:
        hap1 = f"{outdir}/dwgsim_simulated/dwgsim.0.12.fastq",
        hap2 = f"{outdir}/dwgsim_simulated/dwgsim.1.12.fastq",
        fai1 = outdir + "/workflow/input/hap.0.fasta.fai",
        fai2 = outdir + "/workflow/input/hap.1.fasta.fai",
        barcodes = barcodefile
    output:
        collect(outdir + "/lrsim/sim.{hap}.{ext}", hap = [0,1], ext = ["fp", "manifest"]),
        temp(f"{outdir}/lrsim/.status")
    log:
        f"{outdir}/logs/LRSIM.log"
    params:
        lrsim = f"{outdir}/workflow/scripts/LRSIM_harpy.pl",
        proj_dir = f"{outdir}",
        outdist  = config["outer_distance"],
        dist_sd  = config["distance_sd"],
        n_pairs  = config["read_pairs"],
        mol_len  = config["molecule_length"],
        parts    = config["partitions"],
        mols_per = config["molecules_per_partition"],
        static = "-o 1 -d 2 -u 4"
    threads:
        workflow.cores
    conda:
        f"{envdir}/simulations.yaml"
    shell: 
        """
        perl {params.lrsim} -g {input.hap1},{input.hap2} -p {params.proj_dir}/lrsim/sim \\
            -b {input.barcodes} -r {params.proj_dir} -i {params.outdist} \\
            -s {params.dist_sd} -x {params.n_pairs} -f {params.mol_len} \\
            -t {params.parts} -m {params.mols_per} -z {threads} {params.static} 2> {log}
        """

rule sort_manifest:
    input:
        outdir + "/lrsim/sim.{hap}.manifest"
    output:
        outdir + "/lrsim/sim.{hap}.sort.manifest"
    conda:
        f"{envdir}/simulations.yaml"
    shell:
        "msort -kn1 {input} > {output}"

rule extract_reads:
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

rule convert_to_haplotag:
    input:
        fw = outdir + "/10X/sim_hap{hap}_10x_R1_001.fastq.gz",
        rv = outdir + "/10X/sim_hap{hap}_10x_R2_001.fastq.gz",
        barcodes = barcodefile
    output:
        fw = outdir + "/sim_hap{hap}_haplotag.R1.fq.gz",
        rv = outdir + "/sim_hap{hap}_haplotag.R2.fq.gz"
    log:
        outdir + "/logs/10XtoHaplotag/hap{hap}"
    params:
        lambda wc: f"""{outdir}/sim_hap{wc.get("hap")}_haplotag"""
    container:
        None
    shell:
        "10xtoHaplotag.py -f {input.fw} -r {input.rv} -b {input.barcodes} -p {params} > {log}"

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
        summary.append(f"Barcode file: {barcodefile}")
        dwgsim = "Reads were simulated from the provided genomes using:\n"
        dwgsim += f"\tdwgsim -N {params.dwgreadpairs} -e 0.0001,0.0016 -E 0.0001,0.0016 -d {params.dwgouterdist} -s {params.dwgdistsd} -1 135 -2 151 -H -y 0 -S 0 -c 0 -R 0 -r {params.dwgmutationrate} -F 0 -o 1 -m /dev/null GENO PREFIX"
        summary.append(dwgsim)
        lrsim = "LRSIM was started from step 3 (-u 3) with these parameters:\n"
        lrsim += f"\tLRSIM_harpy.pl -g genome1,genome2 -p {params.lrsproj_dir}/lrsim/sim -b BARCODES -r {params.lrsproj_dir} -i {params.lrsoutdist} -s {params.lrsdist_sd} -x {params.lrsn_pairs} -f {params.lrsmol_len} -t {params.lrsparts} -m {params.lrsmols_per} -z THREADS {params.lrsstatic}"
        summary.append(lrsim)
        bxconvert = "10X style barcodes were converted in haplotag BX:Z tags using:\n"
        bxconvert += "\t10xtoHaplotag.py"
        summary.append(bxconvert)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['workflow_call']}"
        summary.append(sm)
        with open(outdir + "/workflow/simulate.reads.summary", "w") as f:
            f.write("\n\n".join(summary))
