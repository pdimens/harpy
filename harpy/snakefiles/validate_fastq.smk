containerized: "docker://pdimens/harpy:latest"

import os
import re
import logging

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake"]["log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

lr_platform = config["linkedreads"]["type"]
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
    params:
        lr_platform
    container:
        None
    shell: 
        "check_fastq {params} {input} > {output}"

rule check_reverse:
    input:
        get_fq2
    output:
        temp("{sample}.R.log")
    params:
        lr_platform
    container:
        None
    shell: 
        "check_fastq {params} {input} > {output}"

rule concat_results:
    input:
        collect("{sample}.{FR}.log", sample = samplenames, FR = ["F","R"])
    output:
        "validate.fastq.tsv"
    container:
        None
    shell:
        """
        echo -e "file\treads\tnoBX\tbadBX\tbadSamSpec\tbxNotLast" > {output}
        cat {input} | sort -k1 >> {output}
        """

rule configure_report:
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
        data = "validate.fastq.tsv",
        qmd = "workflow/report/validate_fastq.qmd"
    output:
        html = "validate.fastq.html",
        qmd = temp("validate.fastq.qmd")
    log:
        "logs/report.log"
    params:
        lr_platform
    conda:
        "envs/report.yaml"
    retries:
        3
    shell:
        """
        PARAM_INFILE=$(realpath {input.data})
        PARAM_PLATFORM="{params}"
        jupyter-nbconvert --execute --output {output.qmd} --to notebook --ExecutePreprocessor.timeout=-1 --ExecutePreprocessor.kernel_name=ir {input.template} 2> {log}
        """

rule workflow_summary:
    default_target: True
    input:
        "validate.fastq.html"
    run:
        summary = ["The harpy validate fastq workflow ran using these parameters:"]
        valids = "Validations were performed with:\n"
        valids += f"\tcheck_fastq {platform} sample.fastq > sample.txt"
        summary.append(valids)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['snakemake']['relative']}"
        summary.append(sm)
        with open("workflow/validate.fastq.summary", "w") as f:
            f.write("\n\n".join(summary))
