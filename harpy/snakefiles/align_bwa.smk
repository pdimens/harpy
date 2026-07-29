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
        bwa_idx = multiext(workflow_geno, ".fai", ".gzi", ".l2b", ".mbw"),
        fai = f"{workflow_geno}.fai",
        gzi = f"{workflow_geno}.gzi" if genome_zip else []
    log:
        f"{workflow_geno}.preprocess.log"
    params:
        genome_zip
    threads:
        workflow.cores
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

            minibwa index -t {threads} {output.geno} 
        }} 2> {log}
        """

rule align:
    input:
        multiext(workflow_geno, ".l2b", ".mbw"),
        ref   = workflow_geno,
        fastq = get_fq
    output:
        bam = temp("bwa/{sample}.bwa.bam"),
        tmp = temp(directory("bwa/{sample}_tmp"))
    log:
        "logs/bwa/{sample}.bwa.log"
    params:
        RG_tag = lambda wc: "-R \"@RG\\tID:" + wc.get("sample") + "\\tSM:" + wc.get("sample") + "\"",
        static = "-x sr -y" if illumina_old else "-x sr",
        extra = extra
    threads:
        12
    resources:
        tmpdir = lambda wc: f"bwa/{wc.sample}_tmp"
    conda:
        "envs/align.yaml"
    container:
        f"docker://pdimens/harpy:align_{VERSION}"
    shell:
        """
        mkdir -p {resources.tmpdir}
        {{
            minibwa map -t {threads} {params} {input.ref} {input.fastq} |
            samtools collate -T {resources.tmpdir} -O -u - 
        }} 2> {log} > {output.bam}
        """

rule sort:
    retries: 3
    input:
        ref = workflow_geno,
        bam = "bwa/{sample}.bwa.bam"
    output:
        bam = temp("sort/{sample}.sort.bam"),
        stats = "reports/data/samtools_stats/{sample}.raw.stats",
        tmp = temp(directory("sort/{sample}_tmp"))
    log:
        "logs/sort/{sample}.sort.log"
    params:
        sortthreads = lambda wc, threads: threads - 1
    threads:
        4
    resources:
        tmpdir = lambda wc: f"sort/{wc.sample}_tmp",
        mem_mb_per_thread = lambda wc, attempt: 3000 // attempt
    shell:
        """
        mkdir -p {resources.tmpdir}
        {{
            samtools fixmate -z on -m -u {input.bam} - |
            samtools sort -@ {params.sortthreads} -M -T {resources.tmpdir} -o {output.bam} -u -l 0 -m {resources.mem_mb_per_thread}M -
            samtools stats -@ {params.sortthreads} -d -x -r {input.ref} {output.bam} > {output.stats} 
        }} 2> {log}
        """

rule mark_duplicates:
    priority: 1
    input:
        fq  = get_fq,
        bam = "sort/{sample}.sort.bam"
    output:
        bam   = "{sample}.bam" if lr_type == "none" or (bx_tag and vx_tag) else temp("markdup/{sample}.bam"),
        stats = "reports/data/markdup/{sample}.markdup",
        tmp = temp(directory("markdup/{sample}_tmp"))
    log:
        "logs/markdup/{sample}.markdup.log"
    params:
        bx_mode = "-S --barcode-tag BX" if not ignore_bx else "-S",
        quality = PARAMETERS.get('min-map-quality', 30),
        unmapped = "-F 4" if not keep_unmapped else "",
        mdthreads = lambda wc, threads: threads - 1
    resources:
        tmpdir = lambda wc: f"markdup/{wc.sample}_tmp"
    threads:
        4
    shell:
        """
        mkdir -p {resources.tmpdir}
        OPT=$(harpy-utils optical-dist-fq {input.fq})
        {{
            samtools view -h -u -q {params.quality} {params.unmapped} {input.bam} |
            samtools markdup -@ {params.mdthreads} -T {resources.tmpdir} {params.bx_mode} -d $OPT -f {output.stats} - {output.bam}
        }} 2> {log}
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

rule depth_stats:
    input:
        "{sample}.bam.bai",
        bam = "{sample}.bam"
    output: 
        "reports/data/coverage/{sample}.regions.bed.gz"
    params:
        f"-b {windowsize}",
        "-n --fast-mode"
    log:
        "logs/depthstats/{sample}.mosdepth.log"
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

rule sample_stats:
    input:
        "{sample}.bam"
    output: 
        temp("{sample}.bam.bai"),
        stats = "reports/data/samtools_stats/{sample}.filtered.stats"
    log:
        "logs/stats/{sample}.stats.log"
    threads:
        2
    shell:
        """
        {{
            samtools index {input}
            samtools stats -@ 1 -x -d {input} > {output.stats}
        }} 2> {log}
        """

rule molecule_coverage:
    input:
        fai = f"{workflow_geno}.fai",
        stats = "reports/data/lrstats/{sample}.lrstats.gz"
    output:
        "reports/data/coverage/{sample}.molcov.gz"
    log:
        "logs/stats/{sample}.molstats.log"
    params:
        windowsize
    shell:
        "harpy-utils molecule-coverage -w {params} {input} 2> {log} | gzip > {output}"

rule molecule_stats:
    input:
        "{sample}.bam"
    output: 
        "reports/data/lrstats/{sample}.lrstats.gz"
    log:
        "logs/molcov/{sample}.molcov.log"
    params:
        molecule_distance
    shell:
        "harpy-utils bx-stats-sam -d {params} {input} 2> {log} | gzip > {output}"

rule alignment_report:
    input:
        collect("reports/data/markdup/{sample}.markdup", sample = samplenames),
        collect("reports/data/samtools_stats/{sample}.{data}.stats", sample = samplenames, data = ["raw", "filtered"]),
        collect("reports/data/coverage/{sample}.regions.bed.gz", sample = samplenames),
        ipynb = f"workflow/samtools_stats.ipynb"
    output:
        tmp = temp("reports/bwa.summary.tmp.ipynb"),
        ipynb = "reports/bwa.summary.ipynb"
    params:
        indir = "-p indir " + os.path.abspath("reports/data")
    log:
        f"logs/reports/bwa.report.log"
    shell:
        """
        export IPYTHONDIR=/tmp/ipython-bwa-stats
        {{
            papermill -k xpython --no-progress-bar --log-level ERROR {input.ipynb} {output.tmp} {params.indir}
            harpy-utils process-notebook {output.tmp} > {output.ipynb}
        }} 2> {log}
        """

rule sample_reports:
    input:
        lrstats = "reports/data/lrstats/{sample}.lrstats.gz",
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
        "logs/reports/{sample}.report.log"
    shell:
        """
        export IPYTHONDIR=/tmp/ipython-{wildcards.sample}.rpt
        {{
            papermill -k xpython --no-progress-bar --log-level ERROR {input.ipynb} {output.tmp} -p platform {params}
            harpy-utils process-notebook {output.tmp} {wildcards.sample} minibwa {params.lr_type} > {output.ipynb}
        }} 2> {log}
        """

rule linked_read_report:
    input:
        collect("reports/data/lrstats/{sample}.lrstats.gz", sample = samplenames),
        ipynb = f"workflow/align_lrstats.ipynb"
    output:
        tmp = temp("reports/linkedreads.summary.tmp.ipynb"),
        ipynb = "reports/linkedreads.summary.ipynb"
    params:
        lr_type = lr_type,
        indir = "-p indir " + os.path.abspath("reports/data/lrstats")
    log:
        f"logs/reports/lrstats.report.log"
    shell:
        """
        {{
            export IPYTHONDIR=/tmp/ipython-bwa.lr
            papermill -k xpython --no-progress-bar --log-level ERROR {input.ipynb} {output.tmp} {params.indir}
            harpy-utils process-notebook {output.tmp} {params.lr_type} > {output.ipynb}
        }} 2> {log}
        """

rule all:
    default_target: True
    input:
        bams = collect("{sample}.bam", sample = samplenames),
        reports = collect("reports/{sample}.ipynb", sample = samplenames) if not skip_reports and not ignore_bx else [],
        align_report = "reports/bwa.summary.ipynb" if (not skip_reports and len(samplenames) > 1) else [],
        bx_report = "reports/linkedreads.summary.ipynb" if (not skip_reports and not ignore_bx and len(samplenames) > 1) else []
