containerized: "docker://pdimens/harpy:latest"

import os
import logging

onstart:
    logger.logger.addHandler(logging.FileHandler(config["snakemake_log"]))
onsuccess:
    os.remove(logger.logfile)
onerror:
    os.remove(logger.logfile)

R1 = config["inputs"]["R1"]
R2 = config["inputs"]["R2"]
I1 = config["inputs"]["I1"]
I2 = config["inputs"]["I2"]
samplefile = config["inputs"]["demultiplex_schema"]
skip_reports = config["skip_reports"]
outdir = config["output_directory"]
envdir = os.getcwd() + "/.harpy_envs"

## the barcode log file ##
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

rule link_R1:
    input:
        R1
    output:
        temp(outdir + "/DATA_R1_001.fastq.gz")
    container:
        None
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

rule barcode_segments:
    output:
        temp(collect(outdir + "/BC_{letter}.txt", letter = ["A","C","B","D"]))
    params:
        outdir
    container:
        None
    shell:
        "haplotag_acbd.py {params}"

rule demux_barcodes:
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
    shell:
        """
        mkdir -p {params.logdir}
        cd {params.outdr}
        demuxGen1 DATA_ demux
        mv demux*BC.log logs
        """

rule demux_samples_R1:
    input:
        f"{outdir}/demux_R1_001.fastq.gz"
    output:
        outdir + "/{sample}.F.fq.gz"
    params:
        c_barcode = lambda wc: samples[wc.get("sample")]
    container:
        None
    shell:
        """
        ( zgrep -A3 "A..{params}B..D" {input} | grep -v "^--$" | gzip -q > {output} ) || touch {output}
        """

use rule demux_samples_R1 as demux_samples_R2 with:
    input:
        f"{outdir}/demux_R2_001.fastq.gz"
    output:
        outdir + "/{sample}.R.fq.gz"

rule fastqc_R1:
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

use rule fastqc_R1 as fastqc_R2 with:
    input:
        outdir + "/{sample}.R.fq.gz"
    output: 
        temp(outdir + "/logs/{sample}_R/fastqc_data.txt")
    params:
        lambda wc: f"{outdir}/logs/" + wc.get("sample") + "_R"

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
    shell:
        "multiqc {params} --filename {output} 2> /dev/null"

rule workflow_summary:
    default_target: True
    input:
        fq = collect(outdir + "/{sample}.{FR}.fq.gz", sample = samplenames, FR = ["F", "R"]),
        reports = outdir + "/reports/demultiplex.QC.html" if not skip_reports else []
    run:
        os.makedirs(f"{outdir}/workflow/", exist_ok= True)
        summary_template = f"""
The harpy demultiplex gen1 workflow ran using these parameters:

Haplotag technology: Generation I

The multiplexed input files:
    - {R1}
    - {R2}
    - {I1}
    - {I2}

Barcodes were moved into the read headers using the command:
    demuxGen1 DATA_ demux

The delimited file associating CXX barcodes with samplenames: {samplefile}

QC checks were performed on demultiplexed FASTQ files using:
    falco -skip-report -skip-summary input.fq.gz

The Snakemake workflow was called via command line:
    {config["workflow_call"]}
"""
        with open(outdir + "/workflow/demux.gen1.summary", "w") as f:
            f.write(summary_template)
