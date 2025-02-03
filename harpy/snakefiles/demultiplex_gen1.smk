containerized: "docker://pdimens/harpy:latest"

import os
import logging

outdir = config["output_directory"]
envdir = os.path.join(os.getcwd(), outdir, "workflow", "envs")
samplefile = config["inputs"]["demultiplex_schema"]
skip_reports = config["reports"]["skip"]

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
            # a casual way to ignore empty lines or lines with !=2 fields
            try:
                smpl, bc = i.split()
                d[smpl] = bc
            except ValueError:
                continue
    return d

samples = barcodedict(samplefile)
samplenames = [i for i in samples.keys()]

rule barcode_segments:
    output:
        temp(collect(outdir + "/workflow/segment_{letter}.bc", letter = ["A","C","B","D"]))
    params:
        f"{outdir}/workflow"
    container:
        None
    shell:
        "haplotag_acbd.py {params}"

rule demultiplex:
    input:
        R1 = config["inputs"]["R1"],
        R2 = config["inputs"]["R2"],
        I1 = config["inputs"]["I1"],
        I2 = config["inputs"]["I2"],
        schema = config["inputs"]["demultiplex_schema"],
        segment_a = f"{outdir}/workflow/segment_A.bc",
        segment_b = f"{outdir}/workflow/segment_B.bc",
        segment_c = f"{outdir}/workflow/segment_C.bc",
        segment_d = f"{outdir}/workflow/segment_D.bc",
    output:
        fw = collect(outdir + "/{sample}.R1.fq.gz", sample = samplenames),
        rv = collect(outdir + "/{sample}.R2.fq.gz", sample = samplenames),
        valid = f"{outdir}/logs/valid_barcodes.log",
        invalid = f"{outdir}/logs/invalid_barcodes.log"
    log:
        f"{outdir}/logs/demultiplex.log"
    params:
        outdir = outdir
    conda:
        f"{envdir}/demultiplex.yaml"
    script:
        "scripts/demultiplex_gen1.py"

rule assess_quality:
    input:
        outdir + "/{sample}.R{FR}.fq.gz"
    output: 
        outdir + "/reports/data/{sample}.R{FR}.fastqc"
    log:
        outdir + "/logs/{sample}.R{FR}.qc.log"
    params:
        f"{outdir}/reports/data"
    threads:
        1
    conda:
        f"{envdir}/qc.yaml"
    shell:
        """
        ( falco --quiet --threads {threads} -skip-report -skip-summary -data-filename {output} {input} ) > {log} 2>&1 ||
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
    params:
        R1 = config["inputs"]["R1"],
        R2 = config["inputs"]["R2"],
        I1 = config["inputs"]["I1"],
        I2 = config["inputs"]["I2"]
    run:
        os.makedirs(f"{outdir}/workflow/", exist_ok= True)
        summary = ["The harpy demultiplex workflow ran using these parameters:"]
        summary.append("Linked Read Barcode Design: Generation I")
        inputs = "The multiplexed input files:\n"
        inputs += f"\tread 1: {params.R1}\n"
        inputs += f"\tread 2: {params.R2}\n"
        inputs += f"\tindex 1: {params.I1}\n"
        inputs += f"\tindex 2: {params.I2}"
        inputs += f"Sample demultiplexing schema: {samplefile}"
        summary.append(inputs)
        demux = "Samples were demultiplexed using:\n"
        demux += "\tworkflow/scripts/demultiplex_gen1.py"
        summary.append(demux)
        qc = "QC checks were performed on demultiplexed FASTQ files using:\n"
        qc += "\tfalco -skip-report -skip-summary -data-filename output input.fq.gz"
        summary.append(qc)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['workflow_call']}"
        summary.append(sm)
        with open(outdir + "/workflow/demux.gen1.summary", "w") as f:
            f.write("\n\n".join(summary))
