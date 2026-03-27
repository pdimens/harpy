import os
import re

localrules: all
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

WORKFLOW    = config.get('Workflow') or {}
PARAMETERS = config.get('Parameters') or {}
REPORTS    = WORKFLOW.get("reports") or {} 
INPUTS     = config['Inputs']
VERSION    = WORKFLOW.get('harpy-version', 'latest')

lr_type           = WORKFLOW.get("linkedreads", {}).get("type", 'none')
bx_tag            = WORKFLOW.get("linkedreads", {}).get("standardized", {}).get("BX", False)
vx_tag            = WORKFLOW.get("linkedreads", {}).get("standardized", {}).get("VX", False)
skip_reports      = REPORTS.get("skip", False)
illumina_old      = PARAMETERS.get("illumina-format-old", False)
molecule_distance = PARAMETERS.get("distance-threshold", 0)
keep_unmapped     = PARAMETERS.get("keep-unmapped", False)
extra 		      = PARAMETERS.get("extra", "") 
windowsize        = PARAMETERS.get("depth-windowsize", 50000)
fqlist            = INPUTS["fastq"]
genomefile 	      = INPUTS["reference"]

ignore_bx     = lr_type == "none"
bn 			  = os.path.basename(genomefile)
workflow_geno = f"workflow/reference/{bn}"
genome_zip    = True if bn.lower().endswith(".gz") else False
geno_idx      = f"{workflow_geno}.gzi" if genome_zip else f"{workflow_geno}.fai"
bn_r          = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
samplenames   = {re.sub(bn_r, "", os.path.basename(i), flags = re.IGNORECASE) for i in fqlist}
d             = dict(zip(samplenames, samplenames))

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

rule optical_dist:
    input:
        get_fq
    output:
        temp("logs/optical/{sample}.opt")
    shell:
        "harpy-utils optical-dist-fq {input} > {output}"

rule align:
    input:
        fastq      = get_fq,
        genome     = workflow_geno,
        genome_idx = multiext(workflow_geno, ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")
    output:
        pipe("samples/{sample}/{sample}.sam")
    log:
        "logs/bwa/{sample}.bwa.log"
    params:
        RG_tag = lambda wc: "-R \"@RG\\tID:" + wc.get("sample") + "\\tSM:" + wc.get("sample") + "\"",
        static = "-C -v 2 -T 10" if illumina_old else "-v 2 -T 10",
        unmapped = "-F 4" if not keep_unmapped else "",
        extra = extra
    threads:
        max(1, min(4, workflow.cores - 2))
    conda:
        "envs/align.yaml"
    container:
        f"docker://pdimens/harpy:align_{VERSION}"
    shell:
        """
        {{
            bwa-mem2 mem -t {threads} {params.RG_tag} {params.static} {params.extra} {input.genome} {input.fastq} |
            samtools view -h -u {params.unmapped}
        }} 2> {log} > {output}
        """     

rule mark_duplicates:
    input:
        sam    = "samples/{sample}/{sample}.sam",
        optical ="logs/optical/{sample}.opt"
    output:
        bam = "{sample}.bam" if lr_type == "none" or (bx_tag and vx_tag) else temp("markdup/{sample}.bam") ,
        stats = "logs/markdup/{sample}.markdup.stats"
    log:
        debug = "logs/markdup/{sample}.markdup.log",
    params: 
        bx_mode = "-S --barcode-tag BX" if not ignore_bx else "-S",
        quality = PARAMETERS.get('min-map-quality', 30),
        opt = lambda wc : open(f"logs/optical/{wc.sample}.opt").read().rstrip(),
        tmprefix = lambda wc: f"samples/{wc.sample}/.{wc.sample}"
    resources:
        mem_mb = 2000
    threads:
        2
    shell:
        """
        {{
            samtools collate -T {params.tmprefix}.collate -O -u {input.sam} - |
            samtools fixmate -z on -m -u - - |
            samtools view -h -u -q {params.quality} |
            samtools sort -T {params.tmprefix}.sort -u -l 0 -m {resources.mem_mb}M - |
            samtools markdup -@ 1 -T {params.tmprefix}.mkdup {params.bx_mode} -d {params.opt} -f {output.stats} - {output.bam}
        }} 2> {log.debug}
        rm -rf {params.tmprefix}*
        """

if lr_type != "none" and not (bx_tag and vx_tag):
    rule standardize:
        input:
            "markdup/{sample}.bam"
        output:
            "{sample}.bam"
        log:
            "logs/{sample}.std.log"
        threads:
            2
        shell:
            "djinn-standardize --threads {threads} {input} > {output} 2> {log}"

rule sample_stats:
    input:
        "{sample}.bam"
    output: 
        temp("{sample}.bam.bai"),
        stats    = temp("reports/data/samtools_stats/{sample}.stats"),
        flagstat = temp("reports/data/samtools_flagstat/{sample}.flagstat"),
        depth    = "reports/data/coverage/{sample}.regions.bed.gz"
    params:
        f"-b {windowsize}",
        "-n --fast-mode"
    log:
        "logs/stats/{sample}.stats.log"
    threads:
        2
    conda:
        "envs/qc.yaml"
    container:
        f"docker://pdimens/harpy:qc_{VERSION}"
    shell:
        """
        {{
            samtools index {input}
            samtools stats -d {input} > {output.stats}
            samtools flagstat {input} > {output.flagstat}
            mosdepth {params} -t 1 reports/data/coverage/{wildcards.sample} {input}
        }} 2> {log}
        rm -f reports/data/coverage/{wildcards.sample}.mosdepth* reports/data/coverage/{wildcards.sample}*.csi
        """

rule molecule_stats:
    input:
        bam = "{sample}.bam",
        fai = f"{workflow_geno}.fai"
    output: 
        stats = "reports/data/bxstats/{sample}.bxstats.gz",
        molcov = "reports/data/coverage/{sample}.molcov.gz"
    log:
        stats = "logs/molcov/{sample}.molcov.log",
        molcov = "logs/stats/{sample}.molstats.log"
    params:
        dist = molecule_distance,
        window = windowsize
    shell:
        """
        harpy-utils bx-stats-sam -d {params.dist} {input.bam} 2> {log.stats} | gzip > {output.stats}
        harpy-utils molecule-coverage -w {params.window} {input.fai} {output.stats} 2> {log.molcov} | gzip > {output.molcov}
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
            papermill -k xpython --no-progress-bar --log-level ERROR {input.ipynb} {output.tmp} -p platform {params}
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
            papermill -k xpython --no-progress-bar --log-level ERROR {input.ipynb} {output.tmp} {params.indir}
            harpy-utils process-notebook {output.tmp} {params.lr_type}
        }} 2> {log} > {output.ipynb}
        """

rule all:
    default_target: True
    input:
        bams = collect("{sample}.bam", sample = samplenames),
        samtools = "reports/bwa.stats.html" if not skip_reports else [],
        reports = collect("reports/{sample}.ipynb", sample = samplenames) if not skip_reports and not ignore_bx else [],
        bx_report = "reports/barcode.summary.ipynb" if (not skip_reports and not ignore_bx and len(samplenames) > 1) else []
