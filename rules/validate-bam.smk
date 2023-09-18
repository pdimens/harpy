import os
import re
import glob

seq_dir = config["seq_directory"]
out_dir_ext = os.path.realpath(config["seq_directory"])
out_dir = "Validate/bam/" + out_dir_ext.replace("/","__") + "/"

bamlist = [os.path.basename(i) for i in glob.iglob(f"{seq_dir}/*") if not os.path.isdir(i) and i.lower().endswith(".bam")]
samplenames = set([os.path.splitext(i)[0] for i in bamlist])  

def get_bam(wildcards):
    # code that returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    lst = [i for i in glob.iglob(seq_dir + "/" + wildcards.sample + "*") if i.lower().endswith(".bam")]
    return lst

def get_bai(wildcards):
    # code that returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    lst = [i for i in glob.iglob(seq_dir + "/" + wildcards.sample + "*") if i.lower().endswith(".bam.bai")]
    return lst

rule validateBam:
    input:
        bam = get_bam,
        bai = get_bai
    output:
        temp(out_dir + "{sample}.log")
    message:
        "Processing: {wildcards.sample}"
    shell: 
        "validateBAM.py {input.bam} > {output}"

rule mergeValidations:
    default_target: True
    input:
        expand(out_dir + "{sample}.log", sample = samplenames)
    output:
        tmp = temp(out_dir + "validations.tmp"),
        final = out_dir + "validations.bam.tsv"
    message:
        "Concatenating results"
    shell:
        """
        cat {input} | sort -k1 > {output.tmp} 
        echo -e "file\tnameMismatch\talignments\tnoBX\tbadBX\n$(cat {output.tmp})" > {output.final}
        """

# TODO add rmd report
#rule createReport:
