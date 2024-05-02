from rich import print as rprint
from rich.panel import Panel
import os
import re
import sys
import glob

seq_dir = config["seq_directory"]
out_dir = config["output_directory"]
envdir      = os.getcwd() + "/.harpy_envs"

bamlist = [os.path.basename(i) for i in glob.iglob(f"{seq_dir}/*") if not os.path.isdir(i) and i.lower().endswith(".bam")]
samplenames = set([os.path.splitext(i)[0] for i in bamlist])  

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

rule index_bam:
    input:
        seq_dir + "/{sample}.bam"
    output:
        seq_dir + "/{sample}.bam.bai"
    message:
        "Indexing {input}"
    shell:
        "samtools index {input}"

rule check_bam:
    input:
        bam = seq_dir + "/{sample}.bam",
        bai = seq_dir + "/{sample}.bam.bai"
    output:
        temp(out_dir + "/{sample}.log")
    conda:
        f"{envdir}/qc.yaml"
    message:
        "Processing: {wildcards.sample}"
    script: 
        "scripts/checkBAM.py"

rule merge_checks:
    input:
        collect(out_dir + "/{sample}.log", sample = samplenames)
    output:
        tmp = temp(out_dir + "/filecheck.tmp"),
        final = out_dir + "/filecheck.bam.tsv"
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
    params:
        seq_dir
    conda:
        f"{envdir}/r.yaml"
    message:
        "Producing report"
    script:
        "report/PreflightBam.Rmd"

rule log_workflow:
    default_target: True
    input:
        out_dir + "/filecheck.bam.html"
    message:
        "Summarizing the workflow: {output}"
    run:
        os.makedirs(f"{out_dir}/workflow/", exist_ok= True)
        with open(out_dir + "/workflow/preflight.bam.summary", "w") as f:
            _ = f.write("The harpy preflight module ran using these parameters:\n\n")
            _ = f.write(f"The directory with sequences: {seq_dir}\n")
            _ = f.write("validations were performed with:\n")
            _ = f.write("    checkBAM.py sample.bam > sample.txt\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")