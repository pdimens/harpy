containerized: "docker://pdimens/harpy:latest"

import os
import re
import logging

onstart:
    logger.logger.addHandler(logging.FileHandler(config["snakemake_log"]))
onsuccess:
    os.remove(logger.logfile)
onerror:
    os.remove(logger.logfile)
wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

fqlist = config["inputs"]
outdir = config["output_directory"]
envdir      = os.path.join(os.getcwd(), outdir, "workflow", "envs")
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
        temp(outdir + "/{sample}.F.log")
    container:
        None
    shell: 
        "check_fastq.py {input} > {output}"

rule check_reverse:
    input:
        get_fq2
    output:
        temp(outdir + "/{sample}.R.log")
    container:
        None
    shell: 
        "check_fastq.py {input} > {output}"

rule concat_results:
    input:
        collect(outdir + "/{sample}.{FR}.log", sample = samplenames, FR = ["F","R"])
    output:
        outdir + "/filecheck.fastq.tsv"
    container:
        None
    shell:
        """
        echo -e "file\treads\tnoBX\tbadBX\tbadSamSpec\tbxNotLast" > {output}
        cat {input} | sort -k1 >> {output}
        """

rule report_config:
    output:
        yaml = f"{outdir}/_quarto.yml",
        scss = f"{outdir}/_harpy.scss"
    params:
        yaml = "https://github.com/pdimens/harpy/raw/refs/heads/main/harpy/reports/_quarto.yml",
        scss = "https://github.com/pdimens/harpy/raw/refs/heads/main/harpy/reports/_harpy.scss"
    run:
        import urllib.request
        with urllib.request.urlopen(params.yaml) as response, open(output.yaml, 'w') as yaml:
            yaml.write(response.read().decode("utf-8"))
        with urllib.request.urlopen(params.scss) as response, open(output.scss, 'w') as scss:
            scss.write(response.read().decode("utf-8"))

rule create_report:
    input:
        data = f"{outdir}/filecheck.fastq.tsv",
        qmd = f"{outdir}/workflow/report/preflight_fastq.qmd",
        f"{outdir}/_quarto.yml",
        f"{outdir}/_harpy.scss"
    output:
        html = f"{outdir}/filecheck.fastq.html",
        qmd = temp(f"{outdir}/filecheck.fastq.qmd")
    log:
        f"{outdir}/logs/report.log"
    conda:
        f"{envdir}/r.yaml"
    shell:
        """
        cp {input.qmd} {output.qmd}
        INFILE=$(realpath {input.data})
        quarto render {output.qmd} --log {log} --quiet -P infile:$INFILE
        """

rule workflow_summary:
    default_target: True
    input:
        outdir + "/filecheck.fastq.html"
    run:
        os.makedirs(f"{outdir}/workflow/", exist_ok= True)
        summary = ["The harpy preflight fastq workflow ran using these parameters:"]
        valids = "Validations were performed with:\n"
        valids += "\tcheck_fastq.py sample.fastq > sample.txt"
        summary.append(valids)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['workflow_call']}"
        summary.append(sm)
        with open(outdir + "/workflow/preflight.fastq.summary", "w") as f:
            f.write("\n\n".join(summary))
