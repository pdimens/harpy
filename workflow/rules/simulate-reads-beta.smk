import os
import re
import sys
import glob
from rich.panel import Panel
from rich import print as rprint

outdir = config["output_directory"]
lrsim_params = f"-p {outdir}/sim"
lrsim_params += " -i " + str(config["outer_distance"])
lrsim_params += " -s " + str(config["distance_sd"])
lrsim_params += " -x " + str(config["read_pairs"])
lrsim_params += " -f " + str(config["molecule_length"])
lrsim_params += " -t " + str(config["partitions"])
lrsim_params += " -m " + str(config["molecules_per_partition"])
gen_hap1 = config["genome_hap1"]
gen_hap2 = config["genome_hap2"]

barcodes = config.get("barcodes", None)
if barcodes:
    barcodefile = barcodes
else:
    barcodefile = f"{outdir}/workflow/input/4M-with-alts-february-2016.txt"

onsuccess:
    print("")
    rprint(
        Panel(
            f"The workflow has finished successfully! Find the results in [bold]{outdir}/[/bold]. If you want to combine both haplotypes of the forward (or reverse) reads together, you can do so with:\n[blue bold]cat reads_hap{{1..2}}.R1.fq.gz > simulations.R1.fq.gz[/blue bold]",
            title = "[bold]harpy simulate reads",
            title_align = "left",
            border_style = "green"
            ),
        file = sys.stderr
    )

onerror:
    print("")
    rprint(
        Panel(
            f"The workflow has terminated due to an error. See the log file below for more details.",
            title = "[bold]harpy simulate reads",
            title_align = "left",
            border_style = "red"
            ),
        file = sys.stderr
    )

rule prepare_directory:
    input:
        hap1 = gen_hap1,
        hap2 = gen_hap2
    output:
        hap1 = f"{outdir}/sim.hap.0.fasta",
        hap2 = f"{outdir}/sim.hap.1.fasta" 
    message:
        "Fooling LRSIM into a false sense of security"
    shell:
        """
        ln -sr {input.hap1} {output.hap1}
        ln -sr {input.hap2} {output.hap2}
        """

rule genome_faidx:
    input:
        outdir + f"{outdir}/sim.hap.{hap}.fasta"
    output: 
        outdir + f"{outdir}/sim.hap.{hap}.fasta.fai"
    log:
        outdir + "/logs/.{hap}.clean.fasta"
    message:
        "Indexing {input}"
    shell:
        "samtools faidx --fai-idx {output} {input} 2> {log}"

if not barcodes:
    rule download_barcodes:
        output:
            barcodefile
        message:
            "Downloading list of standard 10X barcodes"
        run:
            from urllib.request import urlretrieve
            _ = urlretrieve("https://github.com/aquaskyline/LRSIM/blob/master/4M-with-alts-february-2016.txt", output[0])

rule create_molecules_hap1:
    input:
        gen_hap1
    output:
        outdir + "/dwgsim.0.12.fastq"
    params:
        readpairs = config["read_pairs"] * 5000000,
        outerdist = config["outer_distance"],
        distsd = config["distance_sd"]
    conda:
        os.getcwd() + "/.harpy_envs/simulations.yaml"
    message:
        "Creating reads from {input}"
    shell:
        """
        dwgsim -N {params.readpairs} -e 0.0001,0.0016 -E 0.0001,0.0016 -d {params.outerdist} -s {params.distsd} -1 135 -2 151 -H -y 0 -S 0 -c 0 -R 0 -m /dev/null {input} {output}"
        """

use rule create_molecules_hap1 as create_molecules_hap2 with:
    input:
        gen_hap2
    output:
        outdir + "/dwgsim.1.12.fastq"

rule lrsim:
    input:
        hap1 = f"{outdir}/sim.hap.0.fasta",
        hap2 = f"{outdir}/sim.hap.1.fasta",
        fai = expand(outdir + "/sim.hap.{hap}.fasta.fai", hap = [0,1]),
        barcodes = barcodefile
    output:
        expand(outdir + "/sim.{hap}.{ext}", hap = [0,1], ext = ["fp", "manifest"])
    log:
        f"{outdir}/logs/LRSIM.log"
    params:
        lrsim = f"{outdir}/workflow/scripts/LRSIM.pl",
        runoptions = lrsim_params
    threads:
        workflow.cores
    conda:
        os.getcwd() + "/.harpy_envs/simulations.yaml"
    message:
        "Running LRSIM to generate linked reads from\nhaplotype 1: {input.hap1}\nhaplotype 2: {input.hap2}" 
    shell: 
        "{params.lrsim} -g {input.hap1},{input.hap2} {params.runoptions} -z {threads} -o -u 4 > {log}"

rule sort_manifest:
    input:
        outdir + "/sim.{hap}.manifest"
    output:
        outdir + "/sim.{hap}.sort.manifest"
    shell:
        "msort -kn1 {input} > {output}"

rule extract_reads:
    input:
        manifest = outdir + "/sim.{hap}.sort.manifest",
        dwg_hap = outdir + "/dwgsim/hap{hap}.12.fastq"
    output:
        outdir + "/sim_hap{hap}_R1_001.fastq.gz",
        outdir + "/sim_hap{hap}_R1_001.fastq.gz"
    params:
        prefix = lambda wc: outdir + "/sim_hap" + wc.get("hap")
    shell:
        "extractReads {input} {output}"

#TODO adjust python script to be not-snakemake inputs
# add outdir to args
rule convert_haplotag:
    input:
        fw = outdir + "/sim_hap{hap}_R1_001.fastq.gz",
        rv = outdir + "/sim_hap{hap}_R2_001.fastq.gz",
        barcodes = barcodefile
    output:
        fw = outdir + "/sim_hap{hap}_haplotag.R1.fq.gz",
        rv = outdir + "/sim_hap{hap}_haplotag.R2.fq.gz"
    params:
        outdir
    log:
        outdir + "/workflow/10XtoHaplotag_{hap}.txt" 
    message:
        "Converting 10X barcodes to haplotag format"
    shell:
        "10xtoHaplotag.py {input} {outdir}"

rule log_workflow:
    default_target: True
    input:
        expand(outdir + "/hap{hap}_haplotag.R{fw}.fq.gz", hap = [1,2], fw = [1,2])
    output:
        outdir + "/workflow/simulate.reads.workflow.summary"
    message:
        "Creating record of relevant runtime parameters: {output}"
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy simulate reads module ran using these parameters:\n\n")
            _ = f.write(f"Genome haplotype 1: {gen_hap1}\n")
            _ = f.write(f"Genome haplotype 2: {gen_hap2}\n")
            _ = f.write(f"Barcode file: {barcodefile}\n")
            _ = f.write("LRSIM was started from step 3 (-u 3) with these parameters:\n")
            _ = f.write("    " + f"LRSIM.pl -g genome1,genome2 -o -u 3 {lrsim_params}\n")
            _ = f.write("10X style barcodes were converted in haplotag BX:Z tags using:\n")
            _ = f.write("    " + f"10xtoHaplotag.py")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")
