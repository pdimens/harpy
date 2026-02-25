import os
import re

wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

VERSION=4.0
fqlist      = config["Inputs"]["fastq"]
extra 		= config["Parameters"].get("extra", "") 
genomefile 	= config["Inputs"]["reference"]
bn 			= os.path.basename(genomefile)
if bn.lower().endswith(".gz"):
    bn = bn[:-3]
workflow_geno = f"workflow/reference/{bn}"
windowsize  = config["Parameters"]["depth-windowsize"]
molecule_distance = config["Parameters"]["distance-threshold"]
ignore_bx = config["Workflow"]["linkedreads"]["type"] == "none"
is_standardized = config["Workflow"]["linkedreads"]["standardized"]
keep_unmapped = config["Parameters"]["keep-unmapped"]
skip_reports = config["Workflow"]["reports"]["skip"]
plot_contigs = config["Workflow"]["reports"]["plot-contigs"]    
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
        temp("samples/{sample}/{sample}.strobe.sam")
    log:
        "logs/strobealign/{sample}.strobealign.log"
    params: 
        unmapped_strobe = "" if keep_unmapped else "-U",
        unmapped = "" if keep_unmapped else "| samtools view -h -F 4",
        static = "-N 2 -C" if is_standardized else "-N 2",
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
            strobealign {params.static} -t {threads} {params.unmapped_strobe} --rg-id={wildcards.sample} --rg=SM:{wildcards.sample} {params.extra} {input.genome} {input.fastq} {params.unmapped}
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
        cmd = lambda wc: f"samtools collate -O -u samples/{wc.sample}/{wc.sample}.sam" if ignore_bx else f"djinn sam standardize --sam samples/{wc.sample}/{wc.sample}.sam | samtools collate -O -u -",
        bx_mode = "--barcode-tag BX" if not ignore_bx else "",
        quality = config["Parameters"]['min-map-quality']
    resources:
        mem_mb = 2000
    threads:
        2
    shell:
        """
        if grep -q "^[ABCD]" <<< $(samtools head -h 0 -n 1 {input.sam}); then
            OPTICAL_BUFFER=2500
        else
            OPTICAL_BUFFER=100
        fi 
        {{
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
            bx-stats -d {params} {input.bam} | gzip > {output.bxstats}
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
            molecule-coverage -f {input.fai} -w {params} {input.stats} | 
            gzip
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
        "mosdepth {params} -t 1 reports/data/coverage/{wildcards.sample} {input.bam} 2> {log}"

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
        contigs = f"-p contigs {plot_contigs}" if plot_contigs != "default" else "",
        samplename = lambda wc: "-p samplename " + wc.get("sample"),
    log:
        "logs/{sample}.report.log"
    shell:
        """
        {{
            papermill -k python3 --no-progress-bar --log-level ERROR {input.ipynb} {output.tmp} -p platform {params}
            process-notebook {wildcards.sample} strobealign {params.lr_type} {output.tmp}
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
            process-notebook {params.lr_type} {output.tmp}
        }} 2> {log} > {output.ipynb}
        """

rule workflow_summary:
    default_target: True
    input: 
        bams = collect("{sample}.{ext}", sample = samplenames, ext = ["bam","bam.bai"]),
        samtools =  "reports/strobealign.stats.html" if not skip_reports else [] ,
        reports = collect("reports/{sample}.html", sample = samplenames) if not skip_reports and not ignore_bx else [],
        bx_report = "reports/barcode.summary.html" if (not skip_reports and not ignore_bx and len(samplenames) > 1) else []
