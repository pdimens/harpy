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
envdir = os.path.join(os.getcwd(), ".harpy_envs")

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

rule demux_read_1:
    input:
        f"{outdir}/demux_R1_001.fastq.gz"
    output:
        outdir + "/{sample}.R1.fq.gz"
    params:
        c_barcode = lambda wc: samples[wc.get("sample")]
    container:
        None
    shell:
        """
        ( zgrep -A3 "A..{params}B..D" {input} | grep -v "^--$" | gzip -q > {output} ) || touch {output}
        """

use rule demux_read_1 as demux_read_2 with:
    input:
        f"{outdir}/demux_R2_001.fastq.gz"
    output:
        outdir + "/{sample}.R2.fq.gz"

rule quality_assessment:
    input:
        read1 = outdir + "/{sample}.R1.fq.gz",
        read2 = outdir + "/{sample}.R2.fq.gz"
    output:
        temp(outdir + "/reports/data/{sample}.fastp.json")
    log:
        outdir + "/reports/{sample}.html"
    threads:
        2
    conda:
        f"{envdir}/qc.yaml"
    shell:
        "fastp -pQLAG --stdout -w {threads} -i {input.read1} -I {input.read2} -j {output} --html {log} -R \"{wildcards.sample} demultiplex quality report\" > /dev/null"

rule qc_report:
    input:
        collect(outdir + "/reports/data/{sample}.fastp.json", sample = samplenames)
    output:
        outdir + "/reports/demultiplex.QC.html"
    params:
        logdir = outdir + "/reports/data/",
        options = "--no-version-check --force --quiet --no-data-dir",
        title = "--title \"QC for Demultiplexed Samples\"",
        comment = "--comment \"This report aggregates the quality assessments from fastp. The data were NOT filtered, the report is suggesting what it would look like trimmed/filtered with default parameters.\""
    conda:
        f"{envdir}/qc.yaml"
    shell:
        "multiqc {params} --filename {output} 2> /dev/null"

rule workflow_summary:
    default_target: True
    input:
        fq = collect(outdir + "/{sample}.R{FR}.fq.gz", sample = samplenames, FR = [1,2]),
        reports = outdir + "/reports/demultiplex.QC.html" if not skip_reports else []
    run:
        os.makedirs(f"{outdir}/workflow/", exist_ok= True)
        summary = ["The harpy demultiplex workflow ran using these parameters:"]
        summary.append("Haplotag technology: Generation I")
        inputs = "The multiplexed input files:\n"
        inputs += f"\tread 1: {R1}\n"
        inputs += f"\tread 2: {R2}\n"
        inputs += f"\tindex 1: {I1}\n"
        inputs += f"\tindex 2: {I2}"
        summary.append(inputs)
        demux = "Barcodes were moved into the read headers using the command:\n"
        demux += "\tdemuxGen1 DATA_ demux"
        summary.append(demux)
        summary.append(f"The delimited file associating CXX barcodes with samplenames: {samplefile}")
        qc = "QC checks were performed on demultiplexed FASTQ files using:\n"
        qc += "\tfastp -pQLAG --stdout -i R1.fq -I R2.fq > /dev/null"
        summary.append(qc)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['workflow_call']}"
        summary.append(sm)
        with open(outdir + "/workflow/demux.gen1.summary", "w") as f:
            f.write("\n\n".join(summary))
