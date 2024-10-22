containerized: "docker://pdimens/harpy:latest"

import os
import re
import logging
from pathlib import Path

onstart:
    logger.logger.addHandler(logging.FileHandler(config["snakemake_log"]))
onsuccess:
    os.remove(logger.logfile)
onerror:
    os.remove(logger.logfile)
wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

outdir = config["output_directory"]
envdir  = os.path.join(os.getcwd(), ".harpy_envs")
bamlist = config["inputs"]
bamdict = dict(zip(bamlist, bamlist))
samplenames = {Path(i).stem for i in bamlist}

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

rule index_alignments:
    input:
        lambda wc: bamdict[wc.bam]
    output:
        "{bam}.bai"
    container:
        None
    shell:
        "samtools index {input}"

rule check_bam:
    input:
        bam = get_alignments,
        bai = get_align_index
    output:
        temp(outdir + "/{sample}.log")
    container:
        None
    shell: 
        "check_bam.py {input.bam} > {output}"

rule concat_results:
    input:
        collect(outdir + "/{sample}.log", sample = samplenames)
    output:
        tmp = temp(outdir + "/filecheck.tmp"),
        final = outdir + "/filecheck.bam.tsv"
    container:
        None
    shell:
        """
        cat {input} | sort -k1 > {output.tmp} 
        echo -e "file\tnameMismatch\talignments\tnoMI\tnoBX\tbadBX\tbxNotLast\n$(cat {output.tmp})" > {output.final}
        """

rule create_report:
    input:
        outdir + "/filecheck.bam.tsv"
    output:
        outdir + "/filecheck.bam.html"
    log:
        logfile = outdir + "/logs/report.log"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/preflight_bam.Rmd"

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
