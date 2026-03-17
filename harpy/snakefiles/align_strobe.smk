import os
import re

localrules: all
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

WORKFLOW   = config.get('Workflow', {})
PARAMETERS = config.get('Parameters', {})
INPUTS     = config['Inputs']
VERSION    = WORKFLOW.get('harpy-version', 'latest')

lr_type           = WORKFLOW.get("linkedreads", {}).get("type", 'none')
is_standardized   = WORKFLOW.get("linkedreads", {}).get("standardized", False)
skip_reports      = WORKFLOW.get("reports", {}).get("skip", False)
windowsize        = PARAMETERS.get("depth-windowsize", 50000)
molecule_distance = PARAMETERS.get("distance-threshold", 0)
keep_unmapped     = PARAMETERS.get("keep-unmapped", False)
extra 		      = PARAMETERS.get("extra", "") 
fqlist            = INPUTS["fastq"]
genomefile 	      = INPUTS["reference"]

bn 			  = os.path.basename(genomefile)
bn_r          = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
ignore_bx     = lr_type == "none"
bn            = bn[:-3] if bn.lower().endswith(".gz") else bn
workflow_geno = f"workflow/reference/{bn}"
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
        fai = f"{workflow_geno}.fai"
    log:
        f"{workflow_geno}.preprocess.log"
    shell: 
        """
        {{
            seqtk seq {input} > {output.geno}
            samtools faidx --fai-idx {output.fai} {output.geno}
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
        fastq = get_fq,
        genome = workflow_geno
    output:  
        pipe("samples/{sample}/{sample}.strobe.sam")
    log:
        "logs/strobealign/{sample}.strobealign.log"
    params: 
        um_strobe = "" if keep_unmapped else "-U",
        static = "-N 2 -C" if is_standardized else "-N 2",
        RGid = lambda wc: f"--rg-id={wc.get('sample')}",
        RGsm = lambda wc: f"--rg=SM:{wc.get('sample')}",
        extra = extra
    threads:
        min(4, workflow.cores - 2)
    conda:
        "envs/align.yaml"
    container:
        f"docker://pdimens/harpy:align_{VERSION}"
    shell:
        """
        strobealign {params} -t {threads} {input.genome} {input.fastq} 2> {log} > {output}
        """

rule mark_duplicates:
    input:
        sam    = "samples/{sample}/{sample}.strobe.sam",
        genome = workflow_geno,
        faidx  = f"{workflow_geno}.fai",
        optical ="logs/optical/{sample}.opt"
    output:
        "{sample}.bam.bai" if lr_type == "none" or is_standardized else [],
        bam = "{sample}.bam" if lr_type == "none" or is_standardized else temp("markdup/{sample}.bam") 
    log:
        debug = "logs/markdup/{sample}.markdup.log",
        stats = "logs/markdup/{sample}.markdup.stats"
    params: 
        bx_mode = "--barcode-tag BX" if not ignore_bx else "",
        quality = PARAMETERS['min-map-quality'],
        opt = lambda wc : open(f"logs/optical/{wc.sample}.opt").read().rstrip()
    resources:
        mem_mb = 2000
    threads:
        2
    shell:
        """
        {{
            samtools collate -O -u {input.sam} - |
            samtools fixmate -z on -m -u - - |
            samtools view -h -u -q {params.quality} |
            samtools sort -T .{wildcards.sample} -u --reference {input.genome} -l 0 -m {resources.mem_mb}M - |
            samtools markdup -@ 1 -S --write-index {params.bx_mode} -d {params.opt} -f {log.stats} - {output.bam}
        }} 2> {log.debug}
        rm -rf .{wildcards.sample}
        """

if lr_type != "none" or not is_standardized:
    rule standardize:
        input:
            "markdup/{sample}.bam"
        output:
            bai = "{sample}.bam.bai",
            bam = "{sample}.bam"
        log:
            "logs/{sample}.std.log"
        shell:
            """
            {{
                djinn sam standardize {input} > {output.bam}
                samtools index {output.bam}
            }} 2> {log}
            """

rule sample_stats:
    input:
        "{sample}.bam.bai",
        bam = "{sample}.bam"
    output: 
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
            samtools stats -d {input.bam} > {output.stats}
            samtools flagstat {input.bam} > {output.flagstat}
            mosdepth {params} -t 1 reports/data/coverage/{wildcards.sample} {input.bam}
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
        "reports/strobealign.stats.html"
    log:
        "logs/multiqc.log"
    params:
        options = "-n stdout --no-ai  --no-version-check --force --quiet --no-data-dir",
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
        molecule_coverage = "reports/data/coverage/{sample}.molcov.gz",
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
            harpy-utils process-notebook {output.tmp} {wildcards.sample} strobealign {params.lr_type}
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

rule workflow_summary:
    default_target: True
    input: 
        bams = collect("{sample}.{ext}", sample = samplenames, ext = ["bam","bam.bai"]),
        samtools =  "reports/strobealign.stats.html" if not skip_reports else [] ,
        reports = collect("reports/{sample}.ipynb", sample = samplenames) if not skip_reports and not ignore_bx else [],
        bx_report = "reports/barcode.summary.ipynb" if (not skip_reports and not ignore_bx and len(samplenames) > 1) else []
