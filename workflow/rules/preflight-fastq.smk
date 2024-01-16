from rich import print as rprint
from rich.panel import Panel
import os
import re
import sys
import glob

seq_dir = config["seq_directory"]
out_dir = f"{seq_dir}/Preflight/"

flist = [os.path.basename(i) for i in glob.iglob(f"{seq_dir}/*") if not os.path.isdir(i)]
r = re.compile(r".*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
fqlist = list(filter(r.match, flist))
bn_r = r"[\.\_][RF](?:[12])?(?:\_00[1-9])*\.f(?:ast)?q(?:\.gz)?$"
samplenames = set([re.sub(bn_r, "", i, flags = re.IGNORECASE) for i in fqlist])

def get_fq1(wildcards):
    # code that returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    lst = sorted(glob.glob(seq_dir + "/" + wildcards.sample + "*"))
    r = re.compile(r".*[\_\.][FR][1]?(?:\_00[0-9])*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
    fqlist = list(filter(r.match, lst))
    return fqlist[0]

def get_fq2(wildcards):
    # code that returns a list of fastq files for read 2 based on *wildcards.sample*, e.g.
    lst = sorted(glob.glob(seq_dir + "/" + wildcards.sample + "*"))
    r = re.compile(r".*[\_\.][R][2]?(?:\_00[0-9])*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
    fqlist = list(filter(r.match, lst))
    return fqlist

onerror:
    print("")
    rprint(
        Panel(
            f"The workflow has terminated due to an error. See the log file below for more details.",
            title = "[bold]harpy preflight fastq",
            title_align = "left",
            border_style = "red"
            ),
        file = sys.stderr
    )

onsuccess:
    print("")
    rprint(
        Panel(
            f"The workflow has finished successfully! Find the results in [bold]{out_dir}/[/bold]",
            title = "[bold]harpy preflight fastq",
            title_align = "left",
            border_style = "green"
            ),
        file = sys.stderr
    )

rule checkForward:
    input:
        get_fq1
    output:
        temp(out_dir + "{sample}.F.log")
    message:
        "Processing forward reads: {wildcards.sample}"
    shell: 
        "checkFASTQ.py {input} > {output}"

rule checkReverse:
    input:
        get_fq2
    output:
        temp(out_dir + "{sample}.R.log")
    message:
        "Processing reverse reads: {wildcards.sample}"
    shell: 
        "checkFASTQ.py {input} > {output}"

rule mergeChecks:
    input:
        expand(out_dir + "{sample}.{FR}.log", sample = samplenames, FR = ["F","R"])
    output:
        tmp = temp(out_dir + "filecheck.tmp"),
        final = out_dir + "filecheck.fastq.tsv"
    message:
        "Concatenating results"
    shell:
        """
        cat {input} | sort -k1 > {output.tmp}
        echo -e "file\treads\tnoBX\tbadBX\tbadSamSpec\tbxNotLast\n$(cat {output.tmp})" > {output.final}
        """

rule createReport:
    default_target: True
    input:
        out_dir + "filecheck.fastq.tsv"
    output:
        out_dir + "filecheck.fastq.html"
    params:
        seq_dir
    conda:
        os.getcwd() + "/harpyenvs/r-env.yaml"
    message:
        "Producing report"
    script:
        "reportPreflightFastq.Rmd"