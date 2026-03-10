import os
import re
import yaml

wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

VERSION      = config['Workflow']['harpy-version']
fqlist       = config["Inputs"]
skip_reports = config["Workflow"]["reports"]["skip"]
me_seq       = config["Parameters"]["ME-sequence"] 
overlap      = config["Parameters"]["ME-overlap"] 

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
        reports = "reports/preprocess.QA.html" if not skip_reports else []

rule find_ME_seq:
    input:
        get_fq1
    output:
        info = pipe("ME_position/{sample}.info")
    log:
        "logs/{sample}.findME.log"
    params:
        f"-g {me_seq} --overlap {overlap} -e 0.11 --match-read-wildcards --action none -o /dev/null"
    threads:
        2
    conda:
        "envs/qc.yaml"
    container:
        f"docker://pdimens/harpy:qc_{VERSION}"
    shell:
        "cutadapt {params} --info-file {output.info} --cores {threads} {input} 2> {log}"

rule pad_barcodes:
    input:
        info = "ME_position/{sample}.info",
        FQ1 = get_fq1,
        FQ2 = get_fq2
    output:
        temp("stagger/{sample}.stagger.bam")
    log:
        "logs/{sample}.stagger.log"
    shell:
        "gih-stagger {input.FQ1} {input.FQ2} < {input.info} > {output} 2> {log}"

rule extract_barcodes:
    input:
        stagger = "stagger/{sample}.stagger.bam",
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
        pheniqs mux --input {input.stagger} --input {input.stagger} --output {output.bam} --quality -c {input.pheniqs_conf} --report {output.json} 2> {log}
        """

rule count_barcodes:
    input:
        "extract/{sample}.bam"
    output:
        "reports/data/{sample}.rxcount"
    run:
        from collections import Counter
        import pysam
        import re
        barcodes = Counter()
        corrected = 0
        reads = 0
        missed = 0
        invalid = re.compile(r'[ACBD]00')
        with pysam.AlignmentFile(input[0], check_sq=False) as infile:
            for record in infile.fetch(until_eof=True):
                rx = record.get_tag("RX")
                if all(["=" in i for i in rx.split('-')]):
                    missed += 1
                    continue
                qx = record.get_tag("QX")
                corrected += (rx != qx)
                if record.is_read1:
                    reads += 1
                    barcodes.update([rx])
                elif record.is_read2 and (not record.is_paired or _bc not in barcodes):
                    reads += 1
                    barcodes.update([rx])
        with open(output[0], 'w') as fout:
            _unique = sum(not invalid.search(i) for i in barcodes)
            fount.write(f"# total:{reads}|missed:{missed}|corrected barcodes:{corrected}|unique{_unique}\n")
            _ = [fout.write(f"{k}\t{v}\n") for k,v in barcodes.items()]

rule format_barcodes:
    input:
        "workflow/pheniqs.config.json",
        "extract/{sample}.bam"
    output:
        fq1 = "{sample}.R1.fq.gz",
        fq2 = "{sample}.R2.fq.gz"
    log:
        "logs/{sample}.format.BX.log"
    params:
        workflow.cores - 1
    threads:
        workflow.cores
    shell:
        "gih-convert {input} | samtools fastq -@ {params} -N -T VX,BX -1 {output.fq1} -2 {output.fq2} 2> {log}"

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