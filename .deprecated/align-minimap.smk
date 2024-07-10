containerized: "docker://pdimens/harpy:latest"

import os
import re
import sys
from rich.panel import Panel
from rich import print as rprint

outdir      = config["output_directory"]
envdir      = os.getcwd() + "/.harpy_envs"
genomefile 	= config["inputs"]["genome"]
fqlist      = config["inputs"]["fastq"]
extra 		= config.get("extra", "") 
bn 			= os.path.basename(genomefile)
genome_zip  = True if bn.lower().endswith(".gz") else False
bn_idx      = f"{bn}.gzi" if genome_zip else f"{bn}.fai"
skipreports = config["skip_reports"]
windowsize  = config["depth_windowsize"]
molecule_distance = config["molecule_distance"]

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

onerror:
    print("")
    rprint(
        Panel(
            f"The workflow has terminated due to an error. See the log file below for more details.",
            title = "[bold]harpy align minimap",
            title_align = "left",
            border_style = "red"
            ),
        file = sys.stderr
    )

onsuccess:
    print("")
    rprint(
        Panel(
            f"The workflow has finished successfully! Find the results in [bold]{outdir}[/bold]",
            title = "[bold]harpy align minimap",
            title_align = "left",
            border_style = "green"
            ),
        file = sys.stderr
    )

bn_r = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
samplenames = {re.sub(bn_r, "", os.path.basename(i), flags = re.IGNORECASE) for i in fqlist}
d = dict(zip(samplenames, samplenames))

def get_fq(wildcards):
    # returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    r = re.compile(fr".*/({re.escape(wildcards.sample)}){bn_r}", flags = re.IGNORECASE)
    return sorted(list(filter(r.match, fqlist))[:2])

rule genome_setup:
    input:
        genomefile
    output: 
        f"Genome/{bn}"
    container:
        None
    message: 
        "copying {input} to Genome/"
    shell: 
        """
        if (file {input} | grep -q compressed ) ;then
            # is regular gzipped, needs to be BGzipped
            zcat {input} | bgzip -c > {output}
        else
            cp -f {input} {output}
        fi
        """

rule genome_faidx:
    input: 
        f"Genome/{bn}"
    output: 
        fai = f"Genome/{bn}.fai",
        gzi = f"Genome/{bn}.gzi" if genome_zip else []
    log:
        f"Genome/{bn}.faidx.log"
    params:
        genome_zip
    container:
        None
    message:
        "Indexing {input}"
    shell: 
        """
        if [ "{params}" = "True" ]; then
            samtools faidx --gzi-idx {output.gzi} --fai-idx {output.fai} {input} 2> {log}
        else
            samtools faidx --fai-idx {output.fai} {input} 2> {log}
        fi
        """

rule genome_index:
    input: 
        f"Genome/{bn}"
    output: 
        temp(f"Genome/{bn}.mmi")
    log:
        f"Genome/{bn}.mmi.log"
    conda:
        f"{envdir}/align.yaml"
    message:
        "Indexing {input}"
    shell: 
        "minimap2 -x sr -d {output} {input} 2> {log}"

rule align:
    input:
        fastq = get_fq,
        genome   = f"Genome/{bn}.mmi"
    output:  
        pipe(outdir + "/samples/{sample}/{sample}.raw.sam")
    log:
        outdir + "/logs/{sample}.minimap.log"
    params: 
        samps = lambda wc: d[wc.get("sample")],
        extra = extra
    benchmark:
        ".Benchmark/Mapping/minimap/align.{sample}.txt"
    threads:
        min(10, workflow.cores) - 2
    conda:
        f"{envdir}/align.yaml"
    message:
        "Aligning sequences: {wildcards.sample}"
    shell:
        """
        minimap2 -ax sr -t {threads} -y --sam-hit-only -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" {input.genome} {input.fastq} > {output} 2> {log}
        """
 
rule quality_filter:
    input:
        outdir + "/samples/{sample}/{sample}.raw.sam"
    output:
        temp(outdir + "/samples/{sample}/{sample}.mm2.sam")
    params: 
        quality = config["quality"]
    container:
        None
    message:
        "Quality filtering alignments: {wildcards.sample}"
    shell:
        "samtools view -h -F 4 -q {params.quality} {input} > {output}"

rule collate:
    input:
        outdir + "/samples/{sample}/{sample}.mm2.sam"
    output:
        temp(outdir + "/samples/{sample}/{sample}.collate.bam")
    container:
        None
    message:
        "Collating alignments: {wildcards.sample}"
    shell:
        "samtools collate -o {output} {input} 2> /dev/null"

rule fix_mates:
    input:
        outdir + "/samples/{sample}/{sample}.collate.bam"
    output:
        temp(outdir + "/samples/{sample}/{sample}.fixmate.bam")
    container:
        None
    message:
        "Fixing mates in alignments: {wildcards.sample}"
    shell:
        "samtools fixmate -m {input} {output} 2> /dev/null"

rule sort_alignments:
    input:
        sam           = outdir + "/samples/{sample}/{sample}.fixmate.bam",
        genome 		  = f"Genome/{bn}",
        genome_samidx = f"Genome/{bn_idx}"
    output:
        bam = temp(outdir + "/samples/{sample}/{sample}.sort.bam"),
        bai = temp(outdir + "/samples/{sample}/{sample}.sort.bam.bai")
    log:
        outdir + "/logs/{sample}.minimap.sort.log"
    params: 
        quality = config["quality"],
        tmpdir = lambda wc: outdir + "/." + d[wc.sample]
    resources:
        mem_mb = 2000
    container:
        None
    message:
        "Sorting alignments: {wildcards.sample}"
    shell:
        """
        samtools sort -T {params.tmpdir} --reference {input.genome} -O bam -l 0 -m {resources.mem_mb}M --write-index -o {output.bam}##idx##{output.bai} {input.sam} 2> {log}
        rm -rf {params.tmpdir}
        """

rule mark_duplicates:
    input:
        outdir + "/samples/{sample}/{sample}.sort.bam"
    output:
        temp(outdir + "/samples/{sample}/{sample}.markdup.bam")
    log:
        outdir + "/logs/{sample}.markdup.log"
    threads:
        2
    container:
        None
    message:
        "Marking duplicates in alignments alignment: {wildcards.sample}"
    shell:
        "samtools markdup -@ {threads} -S --barcode-tag BX -f {log} {input} {output}  2> /dev/null"

rule index_markdups:
    input:
        outdir + "/samples/{sample}/{sample}.markdup.bam"
    output:
        temp(outdir + "/samples/{sample}/{sample}.markdup.bam.bai")
    container:
        None
    message:
        "Indexing duplicate-marked alignments: {wildcards.sample}"
    shell:
        "samtools index {input}"

rule assign_molecules:
    input:
        bam = outdir + "/samples/{sample}/{sample}.markdup.bam",
        bai = outdir + "/samples/{sample}/{sample}.markdup.bam.bai"
    output:
        bam = outdir + "/{sample}.bam",
        bai = outdir + "/{sample}.bam.bai"
    params:
        molecule_distance
    container:
        None
    message:
        "Assigning barcodes to molecules: {wildcards.sample}"
    shell:
        "assign_mi.py -o {output.bam} -c {params} {input.bam}"

rule alignment_bxstats:
    input:
        bam = outdir + "/{sample}.bam",
        bai = outdir + "/{sample}.bam.bai"
    output: 
        outdir + "/reports/data/bxstats/{sample}.bxstats.gz"
    params:
        sample = lambda wc: d[wc.sample]
    container:
        None
    message:
        "Calculating barcode alignment statistics: {wildcards.sample}"
    shell:
        "bx_stats.py -o {output} {input.bam}"

rule alignment_coverage:
    input: 
        bam = outdir + "/{sample}.bam",
        bai = outdir + "/{sample}.bam.bai"
    output: 
        outdir + "/reports/data/coverage/{sample}.cov.gz"
    params:
        windowsize
    container:
        None
    message:
        "Calculating genomic coverage: {wildcards.sample}"
    shell:
        "samtools depth -a {input.bam} | depth_windows.py {params} | gzip > {output}"

rule alignment_report:
    input:
        outdir + "/reports/data/bxstats/{sample}.bxstats.gz",
        outdir + "/reports/data/coverage/{sample}.cov.gz"
    output:	
        outdir + "/reports/{sample}.html"
    params:
        molecule_distance
    conda:
        f"{envdir}/r.yaml"
    message: 
        "Summarizing barcoded alignments: {wildcards.sample}"
    script:
        "report/AlignStats.Rmd"

rule general_alignment_stats:
    input:
        bam      = outdir + "/{sample}.bam",
        bai      = outdir + "/{sample}.bam.bai"
    output: 
        stats    = temp(outdir + "/reports/data/samtools_stats/{sample}.stats"),
        flagstat = temp(outdir + "/reports/data/samtools_flagstat/{sample}.flagstat")
    container:
        None
    message:
        "Calculating alignment stats: {wildcards.sample}"
    shell:
        """
        samtools stats -d {input.bam} > {output.stats}
        samtools flagstat {input.bam} > {output.flagstat}
        """

rule samtools_reports:
    input: 
        collect(outdir + "/reports/data/samtools_{ext}/{sample}.{ext}", sample = samplenames, ext = ["stats", "flagstat"])
    output: 
        outdir + "/reports/minimap.stats.html"
    params:
        outdir
    conda:
        f"{envdir}/qc.yaml"
    message:
        "Summarizing samtools stats and flagstat"
    shell:
        """
        multiqc {params}/reports/data/samtools_stats {params}/reports/data/samtools_flagstat --no-version-check --force --quiet --title "Basic Alignment Statistics" --comment "This report aggregates samtools stats and samtools flagstats results for all alignments. Samtools stats ignores alignments marked as duplicates." --no-data-dir --filename {output} 2> /dev/null
        """

rule log_workflow:
    default_target: True
    input: 
        bams = collect(outdir + "/{sample}.{ext}", sample = samplenames, ext = ["bam","bam.bai"]),
        samtools =  outdir + "/reports/minimap.stats.html" if not skipreports else [] ,
        bx_reports = collect(outdir + "/reports/{sample}.html", sample = samplenames) if not skipreports else []
    params:
        quality = config["quality"],
        extra   = extra
    message:
        "Summarizing the workflow: {output}"
    run:
        with open(outdir + "/workflow/align.minimap.summary", "w") as f:
            _ = f.write("The harpy align minimap workflow ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write("Sequencing were aligned with Minimap2 using:\n")
            _ = f.write(f"    minimap2 -y {params.extra} --sam-hit-only -R \"@RG\\tID:SAMPLE\\tSM:SAMPLE\" genome.mmi forward_reads reverse_reads |\n")
            _ = f.write(f"    samtools view -h -F 4 -q {params.quality} |\n")
            _ = f.write("Duplicates in the alignments were marked following:\n")
            _ = f.write("    samtools collate \n")
            _ = f.write("    samtools fixmate\n")
            _ = f.write("    samtools sort -m 2000M\n")
            _ = f.write("    samtools markdup -S --barcode-tag BX\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")