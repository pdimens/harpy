import os
import re

wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

VERSION      = config['Workflow']['harpy-version']
fqlist       = config["Inputs"]["fastq"]
molecule_distance = config["Parameters"]["distance-threshold"]
ignore_bx = config["Workflow"]["linkedreads"]["type"] == "none"
is_standardized = config["Workflow"]["linkedreads"]["standardized"]
lr_type = config["Workflow"]["linkedreads"]["type"]
keep_unmapped = config["Parameters"]["keep-unmapped"]
extra 		= config["Parameters"].get("extra", "") 
genomefile 	= config["Inputs"]["reference"]
bn 			= os.path.basename(genomefile)
workflow_geno = f"workflow/reference/{bn}"
genome_zip  = True if bn.lower().endswith(".gz") else False
workflow_geno_idx = f"{workflow_geno}.gzi" if genome_zip else f"{workflow_geno}.fai"
skip_reports = config["Workflow"]["reports"]["skip"]
windowsize  = config["Parameters"]["depth-windowsize"]
bn_r = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
samplenames = {re.sub(bn_r, "", os.path.basename(i), flags = re.IGNORECASE) for i in fqlist}
d = dict(zip(samplenames, samplenames))

def get_fq(wildcards):
    # returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    r = re.compile(fr".*/({re.escape(wildcards.sample)}){bn_r}", flags = re.IGNORECASE)
    return sorted(list(filter(r.match, fqlist))[:2])

rule process_reference:
    input:
        genomefile
    output: 
        geno = workflow_geno,
        bwa_idx = multiext(workflow_geno, ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
        fai = f"{workflow_geno}.fai",
        gzi = f"{workflow_geno}.gzi" if genome_zip else []
    log:
        f"{workflow_geno}.preprocess.log"
    params:
        genome_zip
    conda:
        "envs/align.yaml"
    container:
        f"docker://pdimens/harpy:align_{VERSION}"
    shell: 
        """
        {{
            if (file {input} | grep -q compressed ) ;then
                # is regular gzipped, needs to be BGzipped
                seqtk seq {input} | bgzip -c > {output.geno}
            else
                cp -f {input} {output.geno}
            fi

            if [ "{params}" = "True" ]; then
                samtools faidx --gzi-idx {output.gzi} --fai-idx {output.fai} {output.geno}
            else
                samtools faidx --fai-idx {output.fai} {output.geno}
            fi

            bwa-mem2 index {output.geno} 
        }} 2> {log}
        """

rule align:
    input:
        fastq      = get_fq,
        genome     = workflow_geno,
        genome_idx = multiext(workflow_geno, ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
    output:
        temp("samples/{sample}/{sample}.sam")
    log:
        "logs/bwa/{sample}.bwa.log"
    params:
        RG_tag = lambda wc: "-R \"@RG\\tID:" + wc.get("sample") + "\\tSM:" + wc.get("sample") + "\"",
        static = "-C -v 2 -T 10" if is_standardized else "-v 2 -T 10",
        unmapped = "" if keep_unmapped else "| samtools view -h -F 4",
        extra = extra
    threads:
        min(6, workflow.cores - 1)
    conda:
        "envs/align.yaml"
    container:
        f"docker://pdimens/harpy:align_{VERSION}"
    shell:
        """
        {{
            bwa-mem2 mem -t {threads} {params.RG_tag} {params.static} {params.extra} {input.genome} {input.fastq} {params.unmapped}
        }} 2> {log} > {output}
        """     

rule mark_duplicates:
    input:
        sam    = "samples/{sample}/{sample}.sam",
        genome = workflow_geno,
        faidx  = workflow_geno_idx
    output:
        bai = "{sample}.bam.bai",
        bam = "{sample}.bam"
    log:
        debug = "logs/markdup/{sample}.markdup.log",
        stats = "logs/markdup/{sample}.markdup.stats"
    params:
        cmd = lambda wc: f"samtools collate -O -u samples/{wc.sample}/{wc.sample}.sam" if ignore_bx else f"djinn sam standardize --sam samples/{wc.sample}/{wc.sample}.sam | samtools collate -O -u -",
        bx_mode = "--barcode-tag BX" if not ignore_bx else "",
        quality = config["Parameters"]['min-map-quality']
    resources:
        mem_mb = 2000
    threads:
        4
    shell:
        """
        {{
            OPTICAL_BUFFER=$(harpy-utils optical-dist {input.sam})
            {params.cmd} |
                samtools fixmate -z on -m -u - - |
                samtools view -h -q {params.quality} |
                samtools sort -T .{wildcards.sample} -u --reference {input.genome} -l 0 -m {resources.mem_mb}M - |
                samtools markdup -@ {threads} -S --write-index {params.bx_mode} -d $OPTICAL_BUFFER -f {log.stats} - {output.bam}##idx##{output.bai}
        }} 2> {log.debug}
        rm -rf {wildcards.sample}
        """

rule sample_stats:
    input:
        "{sample}.bam.bai",
        bam = "{sample}.bam"
    output: 
        stats    = temp("reports/data/samtools_stats/{sample}.stats"),
        flagstat = temp("reports/data/samtools_flagstat/{sample}.flagstat"),
        bxstats = "reports/data/bxstats/{sample}.bxstats.gz"
    params:
        molecule_distance
    log:
        "logs/stats/{sample}.stats.log"
    shell:
        """
        {{
            samtools stats -d {input.bam} > {output.stats}
            samtools flagstat {input.bam} > {output.flagstat}
            harpy-utils bx-stats-sam -d {params} {input.bam} | gzip > {output.bxstats}
        }} 2> {log}
        """

rule alignment_coverage:
    input: 
        "{sample}.bam.bai",
        bam = "{sample}.bam"
    output: 
        "reports/data/coverage/{sample}.regions.bed.gz"
    log:
        "logs/stats/{sample}.depth.log"
    params:
        f"-b {windowsize}",
        "-n --fast-mode"
    threads:
        2
    conda:
        "envs/qc.yaml"
    container:
        f"docker://pdimens/harpy:qc_{VERSION}"
    shell:
        """
        mosdepth {params} -t 1 reports/data/coverage/{wildcards.sample} {input.bam} 2> {log}
        rm -f reports/data/coverage/{wildcards.sample}.mosdepth* reports/data/coverage/{wildcards.sample}*.csi
        """

rule molecule_coverage:
    input:
        stats = "reports/data/bxstats/{sample}.bxstats.gz",
        fai = f"{workflow_geno}.fai"
    output:
        "reports/data/coverage/{sample}.molcov.gz"
    log:
        "logs/{sample}.molcov.log"
    params:
        windowsize
    shell:
        """
        {{
            harpy-utils molecule-coverage -w {params} {input.fai} {input.stats} | gzip
        }} > {output} 2> {log}
        """

rule samtools_report:
    input: 
        collect("reports/data/samtools_{ext}/{sample}.{ext}", sample = samplenames, ext = ["stats", "flagstat"])
    output: 
        "reports/bwa.stats.html"
    log:
        "logs/multiqc.log"
    params:
        options = "-n stdout --no-ai --no-version-check --force --quiet --no-data-dir",
        title = "--title \"Basic Alignment Statistics\"",
        comment = "--comment \"This report aggregates samtools stats and samtools flagstats results for all alignments. Samtools stats ignores alignments marked as duplicates.\"",
        outdir = "reports/data/samtools_stats reports/data/samtools_flagstat"
    conda:
        "envs/qc.yaml"
    container:
        f"docker://pdimens/harpy:qc_{VERSION}"
    shell:
        "multiqc {params} > {output} 2> {log}"

rule sample_reports:
    input:
        bxstats = "reports/data/bxstats/{sample}.bxstats.gz",
        coverage = "reports/data/coverage/{sample}.regions.bed.gz",
        molcov = "reports/data/coverage/{sample}.molcov.gz",
        ipynb = f"workflow/align_stats.ipynb"
    output:
        tmp = temp("reports/{sample}.tmp.ipynb"),
        ipynb = "reports/{sample}.ipynb"
    params:
        lr_type = lr_type,
        basedir = "-p basedir " + os.path.abspath("reports/data"),
        mol_dist = f"-p mol_dist {molecule_distance}",
        window_size = f"-p windowsize {windowsize}",
        samplename = lambda wc: "-p samplename " + wc.get("sample"),
    log:
        "logs/{sample}.report.log"
    shell:
        """
        {{
            papermill -k python3 --no-progress-bar --log-level ERROR {input.ipynb} {output.tmp} -p platform {params}
            harpy-utils process-notebook {output.tmp} {wildcards.sample} BWA-MEM2 {params.lr_type}
        }} 2> {log} > {output.ipynb}
        """

rule barcode_report:
    input:
        collect("reports/data/bxstats/{sample}.bxstats.gz", sample = samplenames),
        ipynb = f"workflow/align_bxstats.ipynb"
    output:
        tmp = temp("reports/barcode.summary.tmp.ipynb"),
        ipynb = "reports/barcode.summary.ipynb"
    params:
        lr_type = lr_type,
        indir = "-p indir " + os.path.abspath("reports/data/bxstats")
    log:
        f"logs/reports/bxstats.report.log"
    shell:
        """
        {{
            papermill -k python3 --no-progress-bar --log-level ERROR {input.ipynb} {output.tmp} {params.indir}
            process-notebook {output.tmp} {params.lr_type}
        }} 2> {log} > {output.ipynb}
        """

rule all:
    default_target: True
    input:
        bams = collect("{sample}.{ext}", sample = samplenames, ext = ["bam", "bam.bai"]),
        samtools = "reports/bwa.stats.html" if not skip_reports else [],
        reports = collect("reports/{sample}.ipynb", sample = samplenames) if not skip_reports and not ignore_bx else [],
        bx_report = "reports/barcode.summary.ipynb" if (not skip_reports and not ignore_bx and len(samplenames) > 1) else []
