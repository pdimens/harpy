containerized: "docker://pdimens/harpy:latest"

import os
import re
import logging

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake"]["log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

fqlist      = config["inputs"]["fastq"]
extra 		= config.get("extra", "") 
genomefile 	= config["inputs"]["reference"]
bn 			= os.path.basename(genomefile)
if bn.lower().endswith(".gz"):
    bn = bn[:-3]
workflow_geno = f"workflow/reference/{bn}"
windowsize  = config["depth_windowsize"]
molecule_distance = config["barcodes"]["distance_threshold"]
ignore_bx = config["barcodes"]["ignore"]
is_standardized = config["barcodes"]["standardized"]
keep_unmapped = config["keep_unmapped"]
skip_reports = config["reports"]["skip"]
plot_contigs = config["reports"]["plot_contigs"]    
bn_r = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
samplenames = {re.sub(bn_r, "", os.path.basename(i), flags = re.IGNORECASE) for i in fqlist}
d = dict(zip(samplenames, samplenames))

def get_fq(wildcards):
    # returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    r = re.compile(fr".*/({re.escape(wildcards.sample)}){bn_r}", flags = re.IGNORECASE)
    return sorted(list(filter(r.match, fqlist))[:2])

rule preprocess_reference:
    input:
        genomefile
    output: 
        geno = workflow_geno,
        fai = f"{workflow_geno}.fai"
    log:
        f"{workflow_geno}.preprocess.log"
    container:
        None
    shell: 
        """
        seqtk seq {input} > {output.geno}
        samtools faidx --fai-idx {output.fai} {output.geno} 2> {log}
        """

rule make_depth_intervals:
    input:
        fai = f"{workflow_geno}.fai"
    output:
        bed = "reports/data/coverage/coverage.bed"
    run:
        with open(input.fai, "r") as fai, open(output.bed, "w") as bed:
            for line in fai:
                splitline = line.split()
                contig = splitline[0]
                length = int(splitline[1])
                starts = list(range(0, length, windowsize))
                ends = [i - 1 for i in starts[1:]]
                if not ends or ends[-1] != length:
                    ends.append(length)
                for start,end in zip(starts,ends):
                    bed.write(f"{contig}\t{start}\t{end}\n")

rule align:
    input:
        fastq = get_fq,
        genome   = workflow_geno
    output:  
        pipe("samples/{sample}/{sample}.strobe.sam")
    log:
        "logs/strobealign/{sample}.strobealign.log"
    params: 
        samps = lambda wc: d[wc.get("sample")],
        quality = config["alignment_quality"],
        unmapped_strobe = "" if keep_unmapped else "-U",
        unmapped = "" if keep_unmapped else "-F 4",
        static = "-N 2 -C" if is_standardized else "-N 2",
        extra = extra
    threads:
        min(4, workflow.cores - 1)
    conda:
        "envs/align.yaml"
    shell:
        """
        strobealign {params.static} -t {threads} {params.unmapped_strobe} --rg-id={wildcards.sample} --rg=SM:{wildcards.sample} {params.extra} {input.genome} {input.fastq} 2> {log} |
            samtools view -h {params.unmapped} -q {params.quality} > {output} 
        """

rule standardize_barcodes:
    input:
        "samples/{sample}/{sample}.strobe.sam"
    output:
        temp("samples/{sample}/{sample}.standard.sam")
    container:
        None
    shell:
        "standardize_barcodes_sam.py > {output} < {input}"

rule mark_duplicates:
    input:
        sam    = "samples/{sample}/{sample}.standard.sam",
        genome = workflow_geno,
        faidx  = f"{workflow_geno}.fai"
    output:
        temp("samples/{sample}/{sample}.markdup.bam") if not ignore_bx else temp("markdup/{sample}.markdup.bam")
    log:
        "logs/markdup/{sample}.markdup.log"
    params: 
        tmpdir = lambda wc: "." + d[wc.sample],
        bx_mode = "--barcode-tag BX" if not ignore_bx else ""
    resources:
        mem_mb = 2000
    threads:
        2
    container:
        None
    shell:
        """
        if grep -q "^[ABCD]" <<< $(samtools head -h 0 -n 1 {input.sam}); then
            OPTICAL_BUFFER=2500
        else
            OPTICAL_BUFFER=100
        fi
        samtools collate -O -u {input.sam} |
            samtools fixmate -m -u - - |
            samtools sort -T {params.tmpdir} -u --reference {input.genome} -l 0 -m {resources.mem_mb}M - |
            samtools markdup -@ {threads} -S {params.bx_mode} -d $OPTICAL_BUFFER -f {log} - {output}
        rm -rf {params.tmpdir}
        """

rule assign_molecules:
    priority: 100
    input:
        bam = "samples/{sample}/{sample}.markdup.bam",
    output:
        bam = "{sample}.bam",
        bai = "{sample}.bam.bai"
    params:
        molecule_distance
    container:
        None
    shell:
        """
        assign_mi.py -c {params} {input.bam} > {output.bam}
        samtools index {output.bam}
        """

rule barcode_stats:
    input:
        bam = "{sample}.bam",
        bai = "{sample}.bam.bai"
    output: 
        "reports/data/bxstats/{sample}.bxstats.gz"
    params:
        sample = lambda wc: d[wc.sample]
    container:
        None
    shell:
        "bx_stats.py {input.bam} > {output}"

rule molecule_coverage:
    input:
        stats = "reports/data/bxstats/{sample}.bxstats.gz",
        fai = f"{workflow_geno}.fai"
    output: 
        "reports/data/coverage/{sample}.molcov.gz"
    params:
        windowsize
    container:
        None
    shell:
        "molecule_coverage.py -f {input.fai} -w {params} {input.stats} | gzip > {output}"

rule alignment_coverage:
    input: 
        bam = "{sample}.bam",
        bai = "{sample}.bam.bai",
        bed = "reports/data/coverage/coverage.bed"
    output: 
        "reports/data/coverage/{sample}.cov.gz"
    container:
        None
    shell:
        "samtools bedcov -c {input.bed} {input.bam} | awk '{{ $5 = $5 / ($3 + 1 - $2); print }}' | gzip > {output}"

rule report_config:
    input:
        yaml = "workflow/report/_quarto.yml",
        scss = "workflow/report/_harpy.scss"
    output:
        yaml = temp("reports/_quarto.yml"),
        scss = temp("reports/_harpy.scss")
    run:
        import shutil
        for i,o in zip(input,output):
            shutil.copy(i,o)

rule sample_reports:
    input: 
        "reports/_quarto.yml",
        "reports/_harpy.scss",
        bxstats = "reports/data/bxstats/{sample}.bxstats.gz",
        coverage = "reports/data/coverage/{sample}.cov.gz",
        molecule_coverage = "reports/data/coverage/{sample}.molcov.gz",
        qmd = "workflow/report/align_stats.qmd"
    output:
        report = "reports/{sample}.html",
        qmd = temp("reports/{sample}.qmd")
    params:
        mol_dist = f"-P mol_dist:{molecule_distance}",
        window_size = f"-P windowsize:{windowsize}",
        contigs = f"-P contigs:{plot_contigs}",
        samplename = lambda wc: "-P sample:" + wc.get("sample")
    log:
        "logs/reports/{sample}.alignstats.log"
    conda:
        "envs/r.yaml"
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        BXSTATS=$(realpath {input.bxstats})
        COVFILE=$(realpath {input.coverage})
        MOLCOV=$(realpath {input.molecule_coverage})
        quarto render {output.qmd} --log {log} --quiet -P bxstats:$BXSTATS -P coverage:$COVFILE -P molcov:$MOLCOV {params}
        """

if ignore_bx:
    rule index_bam:
        input:
            bam = "markdup/{sample}.markdup.bam"
        output:
            bam = "{sample}.bam",
            bai = "{sample}.bam.bai"
        container:
            None
        shell:
            """
            mv {input.bam} {output.bam}
            samtools index {output.bam}
            """

rule general_stats:
    input:
        bam      = "{sample}.bam",
        bai      = "{sample}.bam.bai"
    output: 
        stats    = temp("reports/data/samtools_stats/{sample}.stats"),
        flagstat = temp("reports/data/samtools_flagstat/{sample}.flagstat")
    container:
        None
    shell:
        """
        samtools stats -d {input.bam} > {output.stats}
        samtools flagstat {input.bam} > {output.flagstat}
        """

rule samtools_report:
    input: 
        collect("reports/data/samtools_{ext}/{sample}.{ext}", sample = samplenames, ext = ["stats", "flagstat"])
    output: 
        "reports/strobealign.stats.html"
    params:
        outdir = "reports/data/samtools_stats reports/data/samtools_flagstat",
        options = "--no-version-check --force --quiet --no-data-dir",
        title = "--title \"Basic Alignment Statistics\"",
        comment = "--comment \"This report aggregates samtools stats and samtools flagstats results for all alignments. Samtools stats ignores alignments marked as duplicates.\""
    conda:
        "envs/qc.yaml"
    shell:
        "multiqc  {params} --filename {output} 2> /dev/null"

rule barcode_report:
    input: 
        "reports/_quarto.yml",
        "reports/_harpy.scss",
        collect("reports/data/bxstats/{sample}.bxstats.gz", sample = samplenames),
        qmd = "workflow/report/align_bxstats.qmd"
    output:
        report = "reports/barcode.summary.html",
        qmd = temp("reports/barcode.summary.qmd")
    params:
        "reports/data/bxstats/"
    log:
        "logs/reports/bxstats.report.log"
    conda:
        "envs/r.yaml"
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        INPATH=$(realpath {params})
        quarto render {output.qmd} --log {log} --quiet -P indir:$INPATH
        """

rule workflow_summary:
    default_target: True
    input: 
        bams = collect("{sample}.{ext}", sample = samplenames, ext = ["bam","bam.bai"]),
        samtools =  "reports/strobealign.stats.html" if not skip_reports else [] ,
        reports = collect("reports/{sample}.html", sample = samplenames) if not skip_reports and not ignore_bx else [],
        bx_report = "reports/barcode.summary.html" if (not skip_reports and not ignore_bx and len(samplenames) > 1) else []
    params:
        quality = config["alignment_quality"],
        unmapped_strobe = "" if keep_unmapped else "-U",
        unmapped = "" if keep_unmapped else "-F 4",
        bx_mode = "--barcode-tag BX" if not ignore_bx else "",
        static = "-C" if is_standardized else "",
        extra   = extra
    run:
        summary = ["The harpy align strobe workflow ran using these parameters:"]
        summary.append(f"The provided genome: {genomefile}")
        align = "Sequences were aligned with strobealign using:\n"
        align += f"\tstrobealign -U {params.static} --rg-id=SAMPLE --rg=SM:SAMPLE {params.extra} genome reads.F.fq reads.R.fq |\n"
        align += f"\t\tsamtools view -h {params.unmapped} -q {params.quality}"
        summary.append(align)
        standardization = "Barcodes were standardized in the aligments using:\n"
        standardization += "\tstandardize_barcodes_sam.py > {output} < {input}"
        summary.append(standardization)
        duplicates = "Duplicates in the alignments were marked following:\n"
        duplicates += "\tsamtools collate |\n"
        duplicates += "\tsamtools fixmate |\n"
        duplicates += f"\tsamtools sort -T SAMPLE --reference {genomefile} -m 2000M |\n"
        duplicates += f"\tsamtools markdup -S {params.bx_mode} -d 100 (2500 for novaseq)"
        summary.append(duplicates)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['snakemake']['relative']}"
        summary.append(sm)
        with open("workflow/align.strobealign.summary", "w") as f:
            f.write("\n\n".join(summary))