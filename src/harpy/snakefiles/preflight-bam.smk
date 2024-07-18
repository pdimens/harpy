containerized: "docker://pdimens/harpy:latest"

import os
import re
import sys
import glob
import multiprocessing
from pathlib import Path
from rich import print as rprint
from rich.panel import Panel

out_dir = config["output_directory"]
envdir  = os.getcwd() + "/.harpy_envs"
bamlist = config["inputs"]
samplenames = {Path(i).stem for i in bamlist}  

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

onerror:
    print("")
    rprint(
        Panel(
            f"The workflow has terminated due to an error. See the log file below for more details.",
            title = "[bold]harpy preflight bam",
            title_align = "left",
            border_style = "red"
            ),
        file = sys.stderr
    )

onsuccess:
    print("")
    rprint(
        Panel(
            f"The workflow has finished successfully! Find the results in [bold]{out_dir}[/bold]",
            title = "[bold]harpy preflight bam",
            title_align = "left",
            border_style = "green"
            ),
        file = sys.stderr
    )

def get_alignments(wildcards):
    """returns a list with the bam file for the sample based on wildcards.sample"""
    r = re.compile(fr".*/({wildcards.sample})\.(bam|sam)$", flags = re.IGNORECASE)
    aln = list(filter(r.match, bamlist))
    return aln[0]

def get_align_index(wildcards):
    """returns a list with the bai index file for the sample based on wildcards.sample"""
    r = re.compile(fr"(.*/{wildcards.sample})\.(bam|sam)$", flags = re.IGNORECASE)
    aln = list(filter(r.match, bamlist))
    return aln[0] + ".bai"

def sam_index(infile):
    """Use Samtools to index an input file, adding .bai to the end of the name"""
    if not os.path.exists(f"{infile}.bai"):
        subprocess.run(f"samtools index {infile} {infile}.bai".split())

# not the ideal way of doing this, but it works
rule index_alignments:
    input:
        bamlist
    output:
        [f"{i}.bai" for i in bamlist]
    threads:
        workflow.cores
    message:
        "Indexing alignment files"
    run:
        with multiprocessing.Pool(processes=threads) as pool:
            pool.map(sam_index, input)

rule check_bam:
    input:
        bam = get_alignments,
        bai = get_align_index
    output:
        temp(out_dir + "/{sample}.log")
    container:
        None
    message:
        "Processing: {wildcards.sample}"
    shell: 
        "check_bam.py {input.bam} > {output}"

rule merge_checks:
    input:
        collect(out_dir + "/{sample}.log", sample = samplenames)
    output:
        tmp = temp(out_dir + "/filecheck.tmp"),
        final = out_dir + "/filecheck.bam.tsv"
    container:
        None
    message:
        "Concatenating results"
    shell:
        """
        cat {input} | sort -k1 > {output.tmp} 
        echo -e "file\tnameMismatch\talignments\tnoMI\tnoBX\tbadBX\tbxNotLast\n$(cat {output.tmp})" > {output.final}
        """

rule create_report:
    input:
        out_dir + "/filecheck.bam.tsv"
    output:
        out_dir + "/filecheck.bam.html"
    conda:
        f"{envdir}/r.yaml"
    message:
        "Producing report"
    script:
        "report/PreflightBam.Rmd"

rule workflow_summary:
    default_target: True
    input:
        out_dir + "/filecheck.bam.html"
    message:
        "Summarizing the workflow: {output}"
    run:
        os.makedirs(f"{out_dir}/workflow/", exist_ok= True)
        with open(out_dir + "/workflow/preflight.bam.summary", "w") as f:
            _ = f.write("The harpy preflight bam workflow ran using these parameters:\n\n")
            _ = f.write("validations were performed with:\n")
            _ = f.write("    check_bam.py sample.bam > sample.txt\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")