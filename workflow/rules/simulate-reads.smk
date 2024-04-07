import os
import re
import sys
import glob
import gzip
import shutil
from pathlib import Path

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
    shutil.rmtree('./_Inline', ignore_errors=True)
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
    shutil.rmtree('./_Inline', ignore_errors=True)
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


rule genome1_link:
    input:
        gen_hap1
    output: 
        f"{outdir}/workflow/input/hap.0.fasta"
    message: 
        "Decompressing {input}" if gen_hap1.lower().endswith("gz") else "Symlinking {input}"
    run:
        if input[0].lower().endswith("gz"):
            with open(input[0], 'rb') as inf, open(output[0], 'w', encoding='utf8') as outf:
                decom_str = gzip.decompress(inf.read()).decode('utf-8')
                outf.write(decom_str)        
        else:
            if not (Path(output[0]).is_symlink() or Path(output[0]).exists()):
                Path(output[0]).symlink_to(Path(input[0]).absolute()) 

rule genome2_link:
    input:
        gen_hap2
    output: 
        f"{outdir}/workflow/input/hap.1.fasta"
    message: 
        "Decompressing {input}" if gen_hap2.lower().endswith("gz") else "Symlinking {input}"
    run:
        if input[0].lower().endswith("gz"):
            with open(input[0], 'rb') as inf, open(output[0], 'w', encoding='utf8') as outf:
                decom_str = gzip.decompress(inf.read()).decode('utf-8')
                outf.write(decom_str)        
        else:
            if not (Path(output[0]).is_symlink() or Path(output[0]).exists()):
                Path(output[0]).symlink_to(Path(input[0]).absolute()) 

rule genome_faidx:
    input:
        outdir + "/workflow/input/hap.{hap}.fasta"
    output: 
        outdir + "/workflow/input/hap.{hap}.fasta.fai"
    message:
        "Indexing haplotype {wildcards.hap}"
    shell:
        "samtools faidx --fai-idx {output} {input}"

if not barcodes:
    rule download_barcodes:
        output:
            barcodefile
        message:
            "Downloading list of standard 10X barcodes"
        run:
            from urllib.request import urlretrieve
            _ = urlretrieve("https://raw.githubusercontent.com/aquaskyline/LRSIM/master/4M-with-alts-february-2016.txt", output[0])

rule create_molecules_hap1:
    input:
        outdir + "/workflow/input/hap.{hap}.fasta"
    output:
        temp(multiext(outdir + "/dwgsim.{hap}.12", ".bwa.read1.fastq.gz" ,".bwa.read2.fastq.gz", ".mutations.txt", ".mutations.vcf"))
    log:
        outdir + "/logs/dwgsim.hap.{hap}.log"
    params:
        readpairs = config["read_pairs"] * 500000,
        outerdist = config["outer_distance"],
        distsd = config["distance_sd"],
        prefix = lambda wc: outdir + "/dwgsim." + wc.get("hap") + ".12"
    conda:
        os.getcwd() + "/.harpy_envs/simulations.yaml"
    message:
        "Creating reads from {input}"
    shell:
        """
        dwgsim -N {params.readpairs} -e 0.0001,0.0016 -E 0.0001,0.0016 -d {params.outerdist} -s {params.distsd} -1 135 -2 151 -H -y 0 -S 0 -c 0 -R 0 -r 0 -F 0 -o 1 -m /dev/null {input} {params.prefix} 2> {log}
        """

rule interleave_dwgsim_output:
    input:
        expand(outdir + "/dwgsim.{{hap}}.12.bwa.read{rd}.fastq.gz", rd = [1,2]) 
    output:
        outdir + "/dwgsim.{hap}.12.fastq"
    message:
        "Decompressing DWGSIM haplotype {wildcards.hap} output"
    shell:
        "seqtk mergepe {input} > {output}"

rule lrsim:
    default_target: True
    input:
        hap1 = f"{outdir}/dwgsim.0.12.fastq",
        hap2 = f"{outdir}/dwgsim.1.12.fastq",
        fai1 = outdir + "/workflow/input/hap.0.fasta.fai",
        fai2 = outdir + "/workflow/input/hap.1.fasta.fai",
        barcodes = barcodefile
    output:
        expand(outdir + "/sim.{hap}.{ext}", hap = [0,1], ext = ["fp", "manifest"])
    log:
        f"{outdir}/logs/LRSIM.log"
    params:
        lrsim = f"{outdir}/workflow/scripts/LRSIMharpy.pl",
        proj_dir = outdir,
        runoptions = lrsim_params
    threads:
        workflow.cores
    conda:
        os.getcwd() + "/.harpy_envs/simulations.yaml"
    message:
        "Running LRSIM to generate linked reads from\nhaplotype 1: {input.hap1}\nhaplotype 2: {input.hap2}" 
    shell: 
        "perl {params.lrsim} -r {params.proj_dir} -g {input.hap1},{input.hap2} -b {input.barcodes} {params.runoptions} -z {threads} -o -u 4 2> {log}"

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
        outdir + "/sim_hap{hap}_R2_001.fastq.gz"
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
    #default_target: True
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
