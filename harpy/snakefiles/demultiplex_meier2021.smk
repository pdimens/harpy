containerized: "docker://pdimens/harpy:latest"

import os
import logging

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake"]["log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+",
    FR = r"[12]",
    part = r"\d{3}"

schemafile = config["inputs"]["demultiplex_schema"]
skip_reports = config["reports"]["skip"]
qxrx = config["retain"]["qx_rx"]
unknown_samples = config["retain"]["samples"]
unknown_barcodes = config["retain"]["barcodes"]

samplenames = set()
duplicates = False
with open(schemafile, "r") as f:
    for i in f:
        line = i.strip()
        if not i or i.startswith("#"):
            continue
        _sample = line.split()[0] 
        if _sample in samplenames:
            duplicates = True
        samplenames.add(_sample)

if unknown_samples:
    samplenames.add("_unknown_samples")
if unknown_barcodes:
    samplenames.add("_unknown_barcodes")

rule barcode_segments:
    output:
        collect("workflow/segment_{letter}.bc", letter = ["A","C","B","D"])
    container:
        None
    shell:
        "haplotag_acbd workflow"

rule demultiplex:
    input:
        R1 = config["inputs"]["R1"],
        R2 = config["inputs"]["R2"],
        I1 = config["inputs"]["I1"],
        I2 = config["inputs"]["I2"],
        segment_a = "workflow/segment_A.bc",
        segment_b = "workflow/segment_B.bc",
        segment_c = "workflow/segment_C.bc",
        segment_d = "workflow/segment_D.bc",
        schema = schemafile
    output:
        collect("{sample}.R{FR}.fq.gz", sample = samplenames, FR = [1,2]),
        bx_info = "logs/demultiplex.barcodes"
    log:
        "logs/demultiplex.log"
    params:
        outdir = "--samples " + os.getcwd(),
        qxrx = "--rx --qx" if qxrx else "",
        unknown_barcodes = "--undetermined-barcodes _unknown_barcodes" if unknown_barcodes else "",
        unknown_samples = "--undetermined-samples _unknown_samples" if unknown_samples else "",
        duplicate_samples = "--multiple-samples-per-barcode" if duplicates else ""
        #bc_per_segment = "--n-modules {bc_per_segment}"
        #bc_len = "--module-size {bc_len}",
    threads:
        workflow.cores
    conda:
        "envs/demultiplex.yaml"
    shell:
        """
        dmox --i1 {input.I1} --i2 {input.I2} --r1 {input.R1} --r2 {input.R2} \
        --ref-a {input.segment_a} --ref-b {input.segment_b} --ref-c {input.segment_c} \
        --ref-d {input.segment_d} --schema {input.schema} \
        --n-writers {threads} {params} \
        --barcodes-table {output.bx_info} 2> {log}
        """

rule assess_quality:
    input:
        "{sample}.R{FR}.fq.gz"
    output: 
        "reports/data/{sample}.R{FR}.fastqc"
    log:
        "logs/{sample}.R{FR}.qc.log"
    threads:
        1
    conda:
        "envs/qc.yaml"
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

rule configure_report:
    output:
        "workflow/multiqc.yaml"
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
                {"Read the Docs": "https://pdimens.github.io/harpy/workflows/demultiplex/"},
                {"Project Homepage": "https://github.com/pdimens/harpy"}
            ]
        }
        with open(output[0], "w", encoding="utf-8") as yml:
            yaml.dump(configs, yml, default_flow_style= False, sort_keys=False, width=float('inf'))

rule quality_report:
    input:
        fqc = collect("reports/data/{sample}.R{FR}.fastqc", sample = samplenames, FR = [1,2]),
        mqc_yaml = "workflow/multiqc.yaml"
    output:
        "reports/demultiplex.QA.html"
    log:
        "logs/multiqc.log"
    params:
        options = "-n stdout --no-ai --no-version-check --force --quiet --no-data-dir",
        module = " --module fastqc",
        logdir = "reports/data/"
    conda:
        "envs/qc.yaml"
    shell:
        "multiqc --config {input.mqc_yaml} {params} > {output} 2> {log}"

rule all:
    default_target: True
    input:
        fq = collect("{sample}.R{FR}.fq.gz", sample = samplenames, FR = [1,2]),
        barcode_logs = "logs/demultiplex.barcodes",
        reports = "reports/demultiplex.QA.html" if not skip_reports else []
