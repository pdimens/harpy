import os
import re
import yaml

localrules: all
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

WORKFLOW   = config.get('Workflow') or {}
PARAMETERS = config.get('Parameters') or {}
REPORTS    = WORKFLOW.get("reports") or {} 
INPUTS     = config['Inputs']
VERSION    = WORKFLOW.get('harpy-version', 'latest')

fqlist       = INPUTS["fastq"]
skip_reports = REPORTS.get("skip", False)
me_seq       = PARAMETERS.get("ME-sequence", "AGATGTGTATAAGAGACAG")
mismatch     = PARAMETERS.get("ME-mismatch", 2) 
qxrx         = PARAMETERS.get("qx-rx", False)
schemafile   = INPUTS.get("schema", "embedded")

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
            {"Submit an issue": "https://github.com/blinkseq/harpy/issues/new/choose"},
            {"Read the Docs": "https://blinkseq.github.io/harpy/workflows/preprocess/"},
            {"Project Homepage": "https://github.com/blinkseq/harpy"}
        ]
    }
    if not skip_reports:
        with open("workflow/multiqc.yaml", "w", encoding="utf-8") as yml:
            yaml.dump(configs, yml, default_flow_style= False, sort_keys=False, width=float('inf'))

rule all:
    default_target: True
    input: 
        'adapters.fasta',
        collect("{sample}.R{FR}.fq.gz", sample = samplenames, FR = [1,2]),
        collect("reports/data/{sample}.bxcount", sample = samplenames),
        "reports/preprocess.QA.html" if not skip_reports else [],
        "reports/performance.ipynb" if not skip_reports else []

rule write_adapters:
    output:
        'adapters.fasta'
    params:
        me_seq
    run:
        with open(output[0], 'w') as f:
            f.write(">Nextera GIH\nCTGTCTCTTATACACATCT\n")
            #f.write(f">ME Sequence\n{me_seq}\n")

rule demultiplex:
    input:
        FQ1 = get_fq1,
        FQ2 = get_fq2,
        schema = schemafile if schemafile != 'embedded' else []
    output:
        fq1 = "{sample}.R1.fq.gz",
        fq2 = "{sample}.R2.fq.gz",
        stats = "reports/data/{sample}.dmxstats"
    log:
        "logs/{sample}.preproc"
    params:
        f'--me {me_seq}',
        f'--max-mismatch {mismatch}',
        '--qx-rx' if qxrx else '',
        f'-s {schemafile}' if schemafile != 'embedded' else ''
    threads:
        4
    shell:
        "demuxolotl --threads {threads} {params} {wildcards.sample} {input.FQ1} {input.FQ2} > {output.stats} 2> {log}"

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
        collect("reports/data/{sample}.{ext}", sample = samplenames, ext = ['bxcount', 'dmxstats']),
        ipynb = "workflow/preproc_stats.ipynb"
    output:
        "reports/performance.ipynb"
    log:
        "logs/performance.report.log"
    params:
        "-p indir " + os.path.abspath("reports/data/")
    shell:
        """
        export IPYTHONDIR=/tmp/ipython-pre-gih
        harpy-utils run-notebook -k ipython-harpy {params} {input.ipynb} > {output} 2> {log}
        """
