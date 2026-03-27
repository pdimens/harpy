import os
import re
import yaml

localrules: all
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

WORKFLOW   = config.get('Workflow') or {}
PARAMETERS = config.get('Parameters') or {}
INPUTS     = config['Inputs']
VERSION    = WORKFLOW.get('harpy-version', 'latest')

fqlist       = INPUTS
skip_reports = REPORTS.get("skip", False)
me_seq       = PARAMETERS.get("ME-sequence", "AGATGTGTATAAGAGACAG")
mismatch     = PARAMETERS.get("ME-mismatch", 1) 
minlen       = PARAMETERS.get("min-length", 10) 

bn_r = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
samplenames = {re.sub(bn_r, "", os.path.basename(i), flags = re.IGNORECASE) for i in fqlist}

def get_fq1(wildcards):
    # returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    r = re.compile(fr"(.*/{re.escape(wildcards.sample)})([_\.]1|[_\.]F|[_\.]R1(?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$", flags = re.IGNORECASE)
    return list(filter(r.match, fqlist))

def get_fq2(wildcards):
    # returns a list of fastq files for read 2 based on *wildcards.sample*, e.g.
    r = re.compile(fr"(.*/{re.escape(wildcards.sample)})([_\.]2|[_\.]R|[_\.]R2(?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$", flags = re.IGNORECASE)
    return list(filter(r.match, fqlist))

onstart:
    configs = {
        "sp": {"fastqc/data": {"fn" : "*.fastqc"}},
        "table_sample_merge": {
            "R1": ".R1",
            "R2": ".R2"
        },
        "title": "Quality Assessment of Preprocessed Samples",
        "subtitle": "This report aggregates the QA results created by falco",
        "report_comment": "Generated as part of the Harpy preprocess workflow",
        "report_header_info": [
            {"Submit an issue": "https://github.com/pdimens/harpy/issues/new/choose"},
            {"Read the Docs": "https://pdimens.github.io/harpy/workflows/preprocess/"},
            {"Project Homepage": "https://github.com/pdimens/harpy"}
        ]
    }
    if not skip_reports:
        with open("workflow/multiqc.yaml", "w", encoding="utf-8") as yml:
            yaml.dump(configs, yml, default_flow_style= False, sort_keys=False, width=float('inf'))

rule all:
    default_target: True
    input: 
        collect("{sample}.R{FR}.fq.gz", sample = samplenames, FR = [1,2]),
        collect("reports/data/{sample}.bxcount", sample = samplenames),
        "reports/preprocess.QA.html" if not skip_reports else [],
        "reports/performance.ipynb" if not skip_reports else []

rule pad_barcodes:
    input:
        FQ1 = get_fq1,
        FQ2 = get_fq2
    output:
        sam = pipe("stagger/{sample}.stagger.sam"),
        stats = "reports/data/{sample}.MEstats"
    log:
        "logs/{sample}.stagger.log"
    params:
        f'--me {me_seq}',
        f'--max-mismatch {mismatch}',
        f'--min-len {minlen}'
    threads:
        3
    shell:
        "gih-stagger --threads {threads} {params} --stats {output.stats} {input.FQ1} {input.FQ2} > {output.sam} 2> {log}"

rule extract_barcodes:
    input:
        stagger = "stagger/{sample}.stagger.sam",
        pheniqs_conf = "workflow/pheniqs.config.json"
    output:
        bam = temp("extract/{sample}.bam"),
        json = "logs/extract/{sample}.json"
    log:
        "logs/{sample}.pheniqs"
    conda:
        "envs/preprocess.yaml"
    container:
        f"docker://pdimens/harpy:preprocess_{VERSION}"
    shell:
        """
        pheniqs mux --output {output.bam} -s --quality -c {input.pheniqs_conf} --report {output.json} 2> {log} < {input.stagger}
        """

rule format_barcodes:
    input:
        "workflow/pheniqs.config.json",
        "extract/{sample}.bam"
    output:
        fq1 = "{sample}.R1.fq.gz",
        fq2 = "{sample}.R2.fq.gz",
        stats = "reports/data/{sample}.BXstats"
    log:
        "logs/{sample}.convert.log"
    threads:
        2
    shell:
        "gih-convert --threads {threads} {input} {output.fq1} {output.fq2} > {output.stats} 2> {log}"

rule barcode_counts:
    input:
        "{sample}.R1.fq.gz"
    output:
        "reports/data/{sample}.bxcount"
    shell:
        "djinn fastq count {input} > {output}"

rule assess_quality:
    input:
        "{sample}.R{FR}.fq.gz"
    output: 
        "reports/data/{sample}.R{FR}.fastqc"
    log:
        "logs/{sample}.R{FR}.qc.log"
    conda:
        "envs/qc.yaml"
    container:
        f"docker://pdimens/harpy:qc_{VERSION}"
    shell:
        "falco --quiet --threads 1 -skip-report -skip-summary -data-filename {output} {input} > {log} 2>&1"

rule quality_report:
    input:
        fqc = collect("reports/data/{sample}.R{FR}.fastqc", sample = samplenames, FR = [1,2])
    output:
        "reports/preprocess.QA.html"
    log:
        "logs/multiqc.log"
    params:
        options = "-n stdout --no-ai --no-version-check --force --quiet --no-data-dir",
        module = " --module fastqc",
        logdir = "reports/data/"
    conda:
        "envs/qc.yaml"
    container:
        f"docker://pdimens/harpy:qc_{VERSION}"
    shell:
        "multiqc --config workflow/multiqc.yaml {params} > {output} 2> {log}"

rule barcode_report:
    input:
        counts = collect("reports/data/{sample}.bxcount", sample = samplenames),
        stats = collect("reports/data/{sample}.BXstats", sample = samplenames),
        me = collect("reports/data/{sample}.MEstats", sample = samplenames),
        ipynb = "workflow/preproc_stats.ipynb"
    output:
        tmp = temp("reports/performance.tmp.ipynb"),
        ipynb = "reports/performance.ipynb"
    log:
        "logs/performance.report.log"
    params:
        indir = "-p indir " + os.path.abspath("reports/data/")
    shell:
        """
        {{
            papermill -k xpython --no-progress-bar --log-level ERROR {input.ipynb} {output.tmp} {params.indir}
            harpy-utils process-notebook {output.tmp}
        }} 2> {log} > {output.ipynb}
        """
