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
        hap2 = gen_hap2,
    output:
        hap1 = f"{outdir}/sim.hap.0.clean.fasta",
        hap2 = f"{outdir}/sim.hap.1.clean.fasta" 
    message:
        "Fooling LRSIM into a false sense of security"
    shell:
        """
        mkdir -p dwgsim
        ln -s {input.hap1} {output.hap1}
        ln -s {input.hap2} {output.hap2}
        """

rule genome_faidx:
    input:
        outdir + "/sim.hap.{hap}.clean.fasta"
    output: 
        outdir + "/sim.hap.{hap}.clean.fasta.fai"
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
            _ = urlretrieve(
                "https://github.com/aquaskyline/LRSIM/blob/master/4M-with-alts-february-2016.txt",
                output[0]
            )

rule lrsim:
    input:
        hap1 = f"{outdir}/sim.hap.0.clean.fasta",
        hap2 = f"{outdir}/sim.hap.1.clean.fasta",
        fai = expand(outdir + "/sim.hap.{hap}.clean.fasta.fai", hap = [0,1]),
        barcodes = barcodefile
    output:
        expand(outdir + "/sim_S1_L00{hap}_R{fr}_001.fastq.gz", hap = [0,1], fr = [0,1]),
        expand(outdir + "/sim.{hap}.{ext}", hap = [0,1], ext = ["fp", "manifest", "sort.manifest"]),
        temp(expand(outdir + "/sim.dwgsim.{hap}.12.fastq", hap = [0,1]))
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
        "{params.lrsim} -g {input.hap1},{input.hap2} {params.runoptions} -z {threads} -o -u 3 > {log}"

rule convert_haplotag:
    input:
        fw = outdir + "/sim_S1_L00{hap}_R1_001.fastq.gz",
        rv = outdir + "/sim_S1_L00{hap}_R2_001.fastq.gz",
        barcodes = barcodefile
    output:
        fw = "hap{hap}_haplotag.R1.fq.gz",
        rv = "hap{hap}_haplotag.R2.fq.gz"
    log:
        conversions = outdir + "/workflow/10XtoHaplotag_{hap}.txt" 
    message:
        "Converting 10X barcodes to haplotag format"
    script:
        "10xtoHaplotag.py"

rule log_workflow:
    default_target: True
    input:
        expand("hap{hap}_haplotag.R{fw}.fq.gz", hap = [1,2], fw = [1,2])
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
