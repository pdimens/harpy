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

rule align:
    input:
        fastq = get_fq,
        genome   = workflow_geno
    output:  
        temp("samples/{sample}/{sample}.strobe.bam")
    log:
        "logs/strobealign/{sample}.strobealign.log"
    params: 
        um_strobe = "" if keep_unmapped else "-U",
        static = "-N 2 -C" if is_standardized else "-N 2",
        RGid = lambda wc: f"--rg-id={wc.get('sample')}",
        RGsm = lambda wc: f"--rg-SM={wc.get('sample')}"
        extra = extra
    threads:
        min(4, workflow.cores - 1)
    conda:
        "envs/align.yaml"
    container:
        f"docker://pdimens/harpy:align_{VERSION}"
    shell:
        """
        {{
            strobealign {params} -t {threads} {input.genome} {input.fastq} {params.unmapped} |
            samtools view -h -O BAM
        }} 2> {log} > {output} 
        """

rule mark_duplicates:
    input:
        sam    = "samples/{sample}/{sample}.strobe.sam",
        genome = workflow_geno,
        faidx  = f"{workflow_geno}.fai"
    output:
        "{sample}.bam.bai",
        bam = "{sample}.bam"
    log:
        debug = "logs/markdup/{sample}.markdup.log",
        stats = "logs/markdup/{sample}.markdup.stats"
    params: 
        cmd = lambda wc: f"samtools collate -O -u samples/{wc.sample}/{wc.sample}.bam" if ignore_bx else f"djinn sam standardize --sam samples/{wc.sample}/{wc.sample}.bam | samtools collate -O -u -",
        bx_mode = "--barcode-tag BX" if not ignore_bx else "",
        quality = PARAMETERS['min-map-quality']
    resources:
        mem_mb = 2000
    threads:
        2
    shell:
        """
        {{
            OPTICAL_BUFFER=$(harpy-utils optical-dist {input.sam})
            {params.cmd} |
                samtools fixmate -z on -m -u - - |
                samtools view -h -q {params.quality} |
                samtools sort -T .{wildcards.sample} -u --reference {input.genome} -l 0 -m {resources.mem_mb}M - |
                samtools markdup -@ {threads} -S --write-index {params.bx_mode} -d $OPTICAL_BUFFER -f {log.stats} - {output.bam}
        }} 2> {log.debug}
        rm -rf .{wildcards.sample}
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

rule molecule_coverage:
    input:
        stats = "reports/data/bxstats/{sample}.bxstats.gz",
        fai = f"{workflow_geno}.fai"
    output: 
        "reports/data/coverage/{sample}.molcov.gz"
    log:
        "logs/molcov/{sample}.molcov.log"
    params:
        windowsize
    shell:
        """
        {{
            harpy-utils molecule-coverage -w {params} {input.fai} {input.stats} | gzip
        }} > {output} 2> {log}
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
        coverage = "reports/data/coverage/{sample}.cov.gz",
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
            process-notebook {output.tmp} {params.lr_type}
        }} 2> {log} > {output.ipynb}
        """

rule workflow_summary:
    default_target: True
    input: 
        bams = collect("{sample}.{ext}", sample = samplenames, ext = ["bam","bam.bai"]),
        samtools =  "reports/strobealign.stats.html" if not skip_reports else [] ,
        reports = collect("reports/{sample}.html", sample = samplenames) if not skip_reports and not ignore_bx else [],
        bx_report = "reports/barcode.summary.html" if (not skip_reports and not ignore_bx and len(samplenames) > 1) else []
