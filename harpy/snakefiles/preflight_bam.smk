containerized: "docker://pdimens/harpy:latest"

import os
import re
import logging
from pathlib import Path

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake_log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

outdir = config["output_directory"]
envdir  = os.path.join(os.getcwd(), outdir, "workflow", "envs")
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
        temp(outdir + "/{sample}.log")
    container:
        None
    shell: 
        "check_bam.py {input} > {output}"

rule concat_results:
    input:
        collect(outdir + "/{sample}.log", sample = samplenames)
    output:
        outdir + "/filecheck.bam.tsv"
    container:
        None
    shell:
        """
        echo -e "file\talignments\tnameMismatch\tnoMI\tnoBX\tbxNotLast\tbadBX" > {output}
        cat {input} | sort -k1 >> {output}
        """

rule report_config:
    input:
        yaml = f"{outdir}/workflow/report/_quarto.yml",
        scss = f"{outdir}/workflow/report/_harpy.scss"
    output:
        yaml = temp(f"{outdir}/_quarto.yml"),
        scss = temp(f"{outdir}/_harpy.scss")
    run:
        import shutil
        for i,o in zip(input,output):
            shutil.copy(i,o)

rule create_report:
    input:
        f"{outdir}/_quarto.yml",
        f"{outdir}/_harpy.scss",
        data = f"{outdir}/filecheck.bam.tsv",
        qmd = f"{outdir}/workflow/report/preflight_bam.qmd"
    output:
        html = f"{outdir}/filecheck.bam.html",
        qmd = temp(f"{outdir}/filecheck.bam.qmd")
    log:
        f"{outdir}/logs/report.log"
    conda:
        f"{envdir}/r.yaml"
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        INFILE=$(realpath {input.data})
        quarto render {output.qmd} --log {log} --quiet -P infile:$INFILE
        """

rule workflow_summary:
    default_target: True
    input:
        outdir + "/filecheck.bam.html"
    run:
        os.makedirs(f"{outdir}/workflow/", exist_ok= True)
        summary = ["The harpy preflight bam workflow ran using these parameters:"]
        valids = "Validations were performed with:\n"
        valids += "\tcheck_bam.py sample.bam > sample.txt"
        summary.append(valids)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['workflow_call']}"
        summary.append(sm)
        with open(outdir + "/workflow/preflight.bam.summary", "w") as f:
            f.write("\n\n".join(summary))
