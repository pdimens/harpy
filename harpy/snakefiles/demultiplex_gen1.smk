containerized: "docker://pdimens/harpy:latest"

import os
import logging

R1 = config["inputs"]["R1"]
R2 = config["inputs"]["R2"]
I1 = config["inputs"]["I1"]
I2 = config["inputs"]["I2"]
samplefile = config["inputs"]["demultiplex_schema"]
skip_reports = config["skip_reports"]
outdir = config["output_directory"]
envdir = os.path.join(os.getcwd(), ".harpy_envs")

onstart:
    logger.logger.addHandler(logging.FileHandler(config["snakemake_log"]))
    os.makedirs(f"{outdir}/reports/data", exist_ok = True)
    os.makedirs(f"{outdir}/logs/demux", exist_ok = True)
onsuccess:
    os.remove(logger.logfile)
onerror:
    os.remove(logger.logfile)

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

rule demultiplex_barcodes:
    input:
        collect(outdir + "/DATA_{IR}{ext}_001.fastq.gz", IR = ["R","I"], ext = [1,2]),
        collect(outdir + "/BC_{letter}.txt", letter = ["A","C","B","D"])
    output:
        temp(collect(outdir + "/demux_R{ext}_001.fastq.gz", ext = [1,2]))
    params:
        outdir
    container:
        None
    shell:
        """
        cd {params}
        demuxGen1 DATA_ demux
        mv demux*BC.log logs
        """

rule demultiplex_samples:
    input:
        outdir + "/demux_R{FR}_001.fastq.gz"
    output:
        outdir + "/{sample}.R{FR}.fq.gz"
    params:
        c_barcode = lambda wc: samples[wc.get("sample")]
    container:
        None
    shell:
        """
        ( zgrep -A3 "A..{params}B..D" {input} | grep -v "^--$" | gzip -q > {output} ) || touch {output}
        """

rule assess_quality:
    input:
        outdir + "/{sample}.R{FR}.fq.gz"
    output: 
        outdir + "/reports/data/{sample}.R{FR}.fastqc"
    params:
        f"{outdir}/reports/data"
    threads:
        1
    conda:
        f"{envdir}/qc.yaml"
    shell:
        """
        ( falco --quiet --threads {threads} -skip-report -skip-summary -data-filename {output} {input} )> /dev/null ||
cat <<EOF > {output}
##Falco	1.2.4
>>Basic Statistics	fail
#Measure	Value
Filename	{wildcards.sample}.R{wildcards.FR}.fq.gz
File type	Conventional base calls
Encoding	Sanger / Illumina 1.9
Total Sequences	0
Sequences flagged as poor quality	0
Sequence length	0
%GC	0
>>END_MODULE
EOF      
        """

rule report_config:
    output:
        outdir + "/workflow/multiqc.yaml"
    run:
        import yaml
        configs = {
            "sp": {"fastqc/data": {"fn" : "*.fastqc"}},
            "table_sample_merge": {
                "R1": ".R1",
                "R2": ".R2"
            },
            "title": "Quality Assessment of Demultiplexed Samples",
            "subtitle": "This report aggregates the QA results created by falco",
            "report_comment": "Generated as part of the Harpy demultiplex workflow",
            "report_header_info": [
                {"Submit an issue": "https://github.com/pdimens/harpy/issues/new/choose"},
                {"Read the Docs": "https://pdimens.github.io/harpy/"},
                {"Project Homepage": "https://github.com/pdimens/harpy"}
            ]
        }
        with open(output[0], "w", encoding="utf-8") as yml:
            yaml.dump(configs, yml, default_flow_style= False, sort_keys=False, width=float('inf'))

rule qa_report:
    input:
        fqc = collect(outdir + "/reports/data/{sample}.R{FR}.fastqc", sample = samplenames, FR = [1,2]),
        mqc_yaml = outdir + "/workflow/multiqc.yaml"
    output:
        outdir + "/reports/demultiplex.QA.html"
    log:
        f"{outdir}/logs/multiqc.log"
    params:
        options = "--no-version-check --force --quiet --no-data-dir",
        module = " --module fastqc",
        logdir = outdir + "/reports/data/"
    conda:
        f"{envdir}/qc.yaml"
    shell:
        "multiqc --filename {output} --config {input.mqc_yaml} {params} 2> {log}"

rule workflow_summary:
    default_target: True
    input:
        fq = collect(outdir + "/{sample}.R{FR}.fq.gz", sample = samplenames, FR = [1,2]),
        reports = outdir + "/reports/demultiplex.QA.html" if not skip_reports else []
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
