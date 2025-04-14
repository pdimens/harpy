containerized: "docker://pdimens/harpy:latest"

import os
import re
import logging

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake_log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

fqlist = config["inputs"]
bn_r = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
samplenames = {re.sub(bn_r, "", os.path.basename(i), flags = re.IGNORECASE) for i in fqlist}

def get_fq1(wildcards):
    # returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    #samples_FR = [i for i in fqlist if wildcards.sample in i]
    r = re.compile(fr"(.*/{re.escape(wildcards.sample)})([_\.]1|[_\.]F|[_\.]R1(?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$", flags = re.IGNORECASE)
    return list(filter(r.match, fqlist))

def get_fq2(wildcards):
    # returns a list of fastq files for read 2 based on *wildcards.sample*, e.g.
    #samples_FR = [i for i in fqlist if wildcards.sample in i]
    r = re.compile(fr"(.*/{re.escape(wildcards.sample)})([_\.]2|[_\.]R|[_\.]R2(?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$", flags = re.IGNORECASE)
    return list(filter(r.match, fqlist))

rule check_forward:
    input:
        get_fq1
    output:
        temp("{sample}.F.log")
    container:
        None
    shell: 
        "check_fastq.py {input} > {output}"

rule check_reverse:
    input:
        get_fq2
    output:
        temp("{sample}.R.log")
    container:
        None
    shell: 
        "check_fastq.py {input} > {output}"

rule concat_results:
    input:
        collect("{sample}.{FR}.log", sample = samplenames, FR = ["F","R"])
    output:
        "filecheck.fastq.tsv"
    container:
        None
    shell:
        """
        echo -e "file\treads\tnoBX\tbadBX\tbadSamSpec\tbxNotLast" > {output}
        cat {input} | sort -k1 >> {output}
        """

rule report_config:
    input:
        yaml = "workflow/report/_quarto.yml",
        scss = "workflow/report/_harpy.scss"
    output:
        yaml = temp("_quarto.yml"),
        scss = temp("_harpy.scss")
    run:
        import shutil
        for i,o in zip(input,output):
            shutil.copy(i,o)

rule create_report:
    input:
        "_quarto.yml",
        "_harpy.scss",
        data = "filecheck.fastq.tsv",
        qmd = "workflow/report/validate_fastq.qmd"
    output:
        html = "filecheck.fastq.html",
        qmd = temp("filecheck.fastq.qmd")
    log:
        "logs/report.log"
    conda:
        "envs/r.yaml"
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        INFILE=$(realpath {input.data})
        quarto render {output.qmd} --log {log} --quiet -P infile:$INFILE
        """

rule workflow_summary:
    default_target: True
    input:
        "filecheck.fastq.html"
    run:
        summary = ["The harpy validate fastq workflow ran using these parameters:"]
        valids = "Validations were performed with:\n"
        valids += "\tcheck_fastq.py sample.fastq > sample.txt"
        summary.append(valids)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['snakemake_command']}"
        summary.append(sm)
        with open("workflow/validate.fastq.summary", "w") as f:
            f.write("\n\n".join(summary))
