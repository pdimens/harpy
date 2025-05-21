containerized: "docker://pdimens/harpy:latest"

import os
import re
import logging
from pathlib import Path

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake"]["log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

lr_platform = config["platform"]
bamlist = config["inputs"]
bamdict = dict(zip(bamlist, bamlist))
samplenames = {Path(i).stem for i in bamlist}

def get_alignments(wildcards):
    """returns a list with the bam file for the sample based on wildcards.sample"""
    r = re.compile(fr".*/({wildcards.sample})\.(bam|sam)$", flags = re.IGNORECASE)
    aln = list(filter(r.match, bamlist))
    return aln[0]

rule check_bam:
    input:
        get_alignments
    output:
        temp("{sample}.log")
    params:
        lr_platform
    container:
        None
    shell: 
        "check_bam.py {params} {input} > {output}"

rule concat_results:
    input:
        collect("{sample}.log", sample = samplenames)
    output:
        "validate.bam.tsv"
    container:
        None
    shell:
        """
        echo -e "file\talignments\tnameMismatch\tnoMI\tnoBX\tbxNotLast\tbadBX" > {output}
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
        data = "validate.bam.tsv",
        qmd = "workflow/report/validate_bam.qmd"
    output:
        html = "validate.bam.html",
        qmd = temp("validate.bam.qmd")
    params:
        lr_platform
    log:
        "logs/report.log"
    conda:
        "envs/r.yaml"
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        INFILE=$(realpath {input.data})
        quarto render {output.qmd} --log {log} --quiet -P infile:$INFILE -P platform:{params}
        """

rule workflow_summary:
    default_target: True
    input:
        "validate.bam.html"
    run:
        summary = ["The harpy validate bam workflow ran using these parameters:"]
        valids = "Validations were performed with:\n"
        valids += f"\tcheck_bam.py {lr_platform} sample.bam > sample.txt"
        summary.append(valids)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['snakemake']['relative']}"
        summary.append(sm)
        with open("workflow/validate.bam.summary", "w") as f:
            f.write("\n\n".join(summary))
