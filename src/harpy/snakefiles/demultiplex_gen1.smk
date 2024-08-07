containerized: "docker://pdimens/harpy:latest"

import os
import re
import sys
import glob
import logging as pylogging
from datetime import datetime
from rich.panel import Panel
from rich import print as rprint

R1 = config["inputs"]["R1"]
R2 = config["inputs"]["R2"]
I1 = config["inputs"]["I1"]
I2 = config["inputs"]["I2"]
samplefile = config["inputs"]["demultiplex_schema"]
skipreports = config["skip_reports"]
outdir = config["output_directory"]
envdir = os.getcwd() + "/.harpy_envs"

## the log file ##
attempts = glob.glob(f"{outdir}/logs/snakemake/*.snakelog")
if not attempts:
    logfile = f"{outdir}/logs/snakemake/demultiplex_gen1.run1." + datetime.now().strftime("%d_%m_%Y") + ".snakelog"
else:
    increment = sorted([int(i.split(".")[1].replace("run","")) for i in attempts])[-1] + 1
    logfile = f"{outdir}/logs/snakemake/demultiplex_gen1.run{increment}." + datetime.now().strftime("%d_%m_%Y") + ".snakelog"

def barcodedict(smpl):
    d = {}
    with open(smpl, "r") as f:
        for i in f.readlines():
            # a casual way to ignore empty lines or lines with >2 fields
            try:
                smpl, bc = i.split()
                d[smpl] = bc
            except:
                continue
    return d

samples = barcodedict(samplefile)
samplenames = [i for i in samples.keys()]

onstart:
    os.makedirs(f"{outdir}/logs/snakemake", exist_ok = True)
    extra_logfile_handler = pylogging.FileHandler(logfile)
    logger.logger.addHandler(extra_logfile_handler)

onerror:
    print("")
    rprint(
        Panel(
            f"The workflow has terminated due to an error. See the log file for more details:\n[bold]{outdir}/logs/snakemake/{dt_string}.snakelog[/bold]",
            title = "[bold]harpy demultiplex",
            title_align = "left",
            border_style = "red"
            ),
        file = sys.stderr
    )

onsuccess:
    print("")
    rprint(
        Panel(
            f"The workflow has finished successfully! Find the results in [bold]{outdir}[/bold]",
            title = "[bold]harpy demultiplex",
            title_align = "left",
            border_style = "green"
            ),
        file = sys.stderr
    )

rule link_R1:
    input:
        R1
    output:
        temp(outdir + "/DATA_R1_001.fastq.gz")
    container:
        None
    message:
        "Linking {input} to output directory"
    shell:
        "ln -sr {input} {output}"

use rule link_R1 as link_R2 with:
    input: R2
    output: temp(outdir + "/DATA_R2_001.fastq.gz")

use rule link_R1 as link_I1 with:
    input: I1
    output: temp(outdir + "/DATA_I1_001.fastq.gz")

use rule link_R1 as link_I2 with:
    input: I2
    output: temp(outdir + "/DATA_I2_001.fastq.gz")

rule bx_files:
    output:
        temp(collect(outdir + "/BC_{letter}.txt", letter = ["A","C","B","D"]))
    params:
        outdir
    container:
        None
    message:
        "Creating the Gen I barcode files for barcode demultiplexing"
    shell:
        "haplotag_acbd.py {params}"

rule demux_bx:
    input:
        collect(outdir + "/DATA_{IR}{ext}_001.fastq.gz", IR = ["R","I"], ext = [1,2]),
        collect(outdir + "/BC_{letter}.txt", letter = ["A","C","B","D"])
    output:
        temp(collect(outdir + "/demux_R{ext}_001.fastq.gz", ext = [1,2]))
    params:
        outdr = outdir,
        logdir = outdir +"/logs/demux"
    container:
        None
    message:
        "Moving barcodes into read headers"
    shell:
        """
        mkdir -p {params.logdir}
        cd {params.outdr}
        demuxGen1 DATA_ demux
        mv demux*BC.log logs
        """

rule split_samples_fw:
    input:
        f"{outdir}/demux_R1_001.fastq.gz"
    output:
        outdir + "/{sample}.F.fq.gz"
    params:
        c_barcode = lambda wc: samples[wc.get("sample")]
    container:
        None
    message:
        "Extracting forward reads:\n sample: {wildcards.sample}\n barcode: {params}"
    shell:
        """
        ( zgrep -A3 "A..{params}B..D" {input} | grep -v "^--$" | gzip -q > {output} ) || touch {output}
        """

use rule split_samples_fw as split_samples_rv with:
    input:
        f"{outdir}/demux_R2_001.fastq.gz"
    output:
        outdir + "/{sample}.R.fq.gz"
    message:
        "Extracting reverse reads:\n sample: {wildcards.sample}\n barcode: {params}"

rule fastqc_F:
    input:
        outdir + "/{sample}.F.fq.gz"
    output: 
        temp(outdir + "/logs/{sample}_F/fastqc_data.txt")
    params:
        lambda wc: f"{outdir}/logs/" + wc.get("sample") + "_F"
    threads:
        1
    conda:
        f"{envdir}/qc.yaml"
    message:
        "Performing quality assessment: {wildcards.sample}.F.fq.gz"
    shell:
        """
        mkdir -p {params}
        if [ -z $(gzip -cd {input} | head -c1) ]; then
            echo "##Falco	1.2.1" > {output}
            echo ">>Basic Statistics	fail" >> {output}
            echo "#Measure	Value" >> {output}
            echo "Filename	{wildcards.sample}.F.fq.gz" >> {output}
            echo "File type	Conventional base calls" >> {output}
            echo "Encoding	Sanger / Illumina 1.9" >> {output}
            echo "Total Sequences	0" >> {output}
            echo "Sequences flagged as poor quality	0" >> {output}
            echo "Sequence length	0" >> {output}
            echo "%GC	0" >> {output}
            echo ">>END_MODULE" >> {output}
        else
            falco -q --threads {threads} -skip-report -skip-summary -o {params} {input}
        fi
        """

use rule fastqc_F as fastqc_R with:
    input:
        outdir + "/{sample}.R.fq.gz"
    output: 
        temp(outdir + "/logs/{sample}_R/fastqc_data.txt")
    params:
        lambda wc: f"{outdir}/logs/" + wc.get("sample") + "_R"
    message:
        "Performing quality assessment: {wildcards.sample}.R.fq.gz"

rule qc_report:
    input:
        collect(outdir + "/logs/{sample}_{FR}/fastqc_data.txt", sample = samplenames, FR = ["F","R"])
    output:
        outdir + "/reports/demultiplex.QC.html"
    params:
        logdir = outdir + "/logs/",
        options = "--no-version-check --force --quiet --no-data-dir",
        title = "--title \"QC for Demultiplexed Samples\"",
        comment = "--comment \"This report aggregates the QC results created by falco.\""
    conda:
        f"{envdir}/qc.yaml"
    message:
        "Creating final demultiplexing QC report"
    shell:
        "multiqc {params} --filename {output} 2> /dev/null"

rule workflow_summary:
    default_target: True
    input:
        fq = collect(outdir + "/{sample}.{FR}.fq.gz", sample = samplenames, FR = ["F", "R"]),
        reports = outdir + "/reports/demultiplex.QC.html" if not skipreports else []
    message:
        "Summarizing the workflow: {output}"
    run:
        os.makedirs(f"{outdir}/workflow/", exist_ok= True)
        with open(outdir + "/workflow/demux.gen1.summary", "w") as f:
            _ = f.write("The harpy demultiplex gen1 workflow ran using these parameters:\n\n")
            _ = f.write("Haplotag technology: Generation I\n")
            _ = f.write(f"The multiplexed input files:\n    -")
            _ = f.write("\n    -".join([R1,R2,I1,I2]) + "\n")
            _ = f.write("Barcodes were moved into the read headers using the command:\n")
            _ = f.write(f"    demuxGen1 DATA_ demux\n")
            _ = f.write(f"The delimited file associating CXX barcodes with samplenames:\n    {samplefile}\n")
            _ = f.write(f"QC checks were performed on demultiplexed FASTQ files using:\n")
            _ = f.write(f"    falco -skip-report -skip-summary input.fq.gz\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")
