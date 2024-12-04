import os
import re
import glob

seq_dir = config["seq_directory"]
out_dir = config["output_directory"]

bamlist = [os.path.basename(i) for i in glob.iglob(f"{seq_dir}/*") if not os.path.isdir(i) and i.lower().endswith(".bam")]
samplenames = set([os.path.splitext(i)[0] for i in bamlist])  

#def get_bam(wildcards):
#    # returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
#    lst = [i for i in glob.iglob(seq_dir + "/" + wildcards.sample + "*") if i.lower().endswith(".bam")]
#    return lst
#
#def get_bai(wildcards):
#    # returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
#    lst = [i for i in glob.iglob(seq_dir + "/" + wildcards.sample + "*") if i.lower().endswith(".bam.bai")]
#    return lst

rule indexBam:
    input:
        seq_dir + "/{sample}.bam"
    output:
        seq_dir + "/{sample}.bam.bai"
    message:
        "Indexing {input}"
    shell:
        "samtools index {input}"

rule checkBam:
    input:
        bam = seq_dir + "/{sample}.bam",
        bai = seq_dir + "/{sample}.bam.bai"
    output:
        temp(out_dir + "{sample}.log")
    conda:
        os.getcwd() + "/.harpy_envs/qc.yaml"
    message:
        "Processing: {wildcards.sample}"
    script: 
        "scripts/checkBAM.py"

rule mergeChecks:
    input:
        expand(out_dir + "{sample}.log", sample = samplenames)
    output:
        tmp = temp(out_dir + "filecheck.tmp"),
        final = out_dir + "filecheck.bam.tsv"
    message:
        "Concatenating results"
    shell:
        """
        cat {input} | sort -k1 > {output.tmp} 
        echo -e "file\tnameMismatch\talignments\tnoMI\tnoBX\tbadBX\tbxNotLast\n$(cat {output.tmp})" > {output.final}
        """

rule createReport:
    default_target: True
    input:
        out_dir + "filecheck.bam.tsv"
    output:
        out_dir + "filecheck.bam.html"
    params:
        seq_dir
    message:
        "Producing report"
    script:
        "reportPreflightBam.Rmd"
