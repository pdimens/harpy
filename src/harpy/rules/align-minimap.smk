import os
import re
import glob
import sys
from rich.panel import Panel
from rich import print as rprint

outdir      = config["output_directory"]
seq_dir		= config["seq_directory"]
genomefile 	= config["genomefile"]
samplenames = config["samplenames"]
molecule_distance = config["molecule_distance"]
extra 		= config.get("extra", "") 
bn 			= os.path.basename(genomefile)
genome_zip  = True if bn.lower().endswith(".gz") else False
bn_idx      = f"{bn}.gzi" if genome_zip else f"{bn}.fai"
skipreports = config["skipreports"]

d = dict(zip(samplenames, samplenames))

def get_fq1(wildcards):
    # returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    # the list is just the single reverse file
    lst = glob.glob(seq_dir + "/" + wildcards.sample + "*")
    return lst[0]

def get_fq2(wildcards):
    # returns a list of fastq files for read 2 based on *wildcards.sample*, e.g.
    # the list is just the single reverse file
    lst = sorted(glob.glob(seq_dir + "/" + wildcards.sample + "*"))
    return lst[1]

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

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

rule genome_link:
    input:
        genomefile
    output: 
        f"Genome/{bn}"
    message: 
        "Symlinking {input}"
    shell: 
        """
        if (file {input} | grep -q compressed ) ;then
            # is regular gzipped, needs to be BGzipped
            zcat {input} | bgzip -c > {output}
        elif (file {input} | grep -q BGZF ); then
            # is bgzipped, just linked
            ln -sr {input} {output}
        else
            # isn't compressed, just linked
            ln -sr {input} {output}
        fi
        """

if genome_zip:
    rule genome_compressed_faidx:
        input: 
            f"Genome/{bn}"
        output: 
            gzi = f"Genome/{bn}.gzi",
            fai = f"Genome/{bn}.fai"
        log:
            f"Genome/{bn}.faidx.gzi.log"
        message:
            "Indexing {input}"
        shell: 
            "samtools faidx --gzi-idx {output.gzi} --fai-idx {output.fai} {input} 2> {log}"
else:
    rule genome_faidx:
        input: 
            f"Genome/{bn}"
        output: 
            f"Genome/{bn}.fai"
        log:
            f"Genome/{bn}.faidx.log"
        message:
            "Indexing {input}"
        shell:
            "samtools faidx --fai-idx {output} {input} 2> {log}"

rule genome_make_windows:
    input:
        f"Genome/{bn}"
    output: 
        f"Genome/{bn}.bed"
    message: 
        "Creating BED intervals from {input}"
    shell: 
        "makewindows.py -i {input} -w 10000 -o {output}"

rule genome_index:
    input: 
        f"Genome/{bn}"
    output: 
        temp(f"Genome/{bn}.mmi")
    log:
        f"Genome/{bn}.mmi.log"
    conda:
        os.getcwd() + "/.harpy_envs/align.yaml"
    message:
        "Indexing {input}"
    shell: 
        "minimap2 -x sr -d {output} {input} 2> {log}"

rule align:
    input:
        forward_reads = get_fq1,
        reverse_reads = get_fq2,
        genome 		  = f"Genome/{bn}.mmi"
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
        os.getcwd() + "/.harpy_envs/align.yaml"
    message:
        "Aligning sequences: {wildcards.sample}"
    shell:
        """
        minimap2 -ax sr  -t {threads} -y --sam-hit-only -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" {input.genome} {input.forward_reads} {input.reverse_reads} > {output} 2> {log}
        """
 
rule quality_filter:
    input:
        outdir + "/samples/{sample}/{sample}.raw.sam"
    output:
        temp(outdir + "/samples/{sample}/{sample}.mm2.sam")
    params: 
        quality = config["quality"]
    message:
        "Quality filtering alignments: {wildcards.sample}"
    shell:
        "samtools view -h -F 4 -q {params.quality} {input} > {output}"

rule collate:
    input:
        outdir + "/samples/{sample}/{sample}.mm2.sam"
    output:
        temp(outdir + "/samples/{sample}/{sample}.collate.bam")
    message:
        "Collating alignments: {wildcards.sample}"
    shell:
        "samtools collate -o {output} {input} 2> /dev/null"

rule fix_mates:
    input:
        outdir + "/samples/{sample}/{sample}.collate.bam"
    output:
        temp(outdir + "/samples/{sample}/{sample}.fixmate.bam")
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
    message:
        "Sorting alignments: {wildcards.sample}"
    shell:
        """
        samtools sort -T {params.tmpdir} --reference {input.genome} -O bam -l 0 -m 4G --write-index -o {output.bam}##idx##{output.bai} {input.sam} 2> {log}
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
    message:
        "Marking duplicates in alignments alignment: {wildcards.sample}"
    shell:
        "samtools markdup -@ {threads} -S --barcode-tag BX -f {log} {input} {output}  2> /dev/null"

rule index_markdups:
    input:
        outdir + "/samples/{sample}/{sample}.markdup.bam"
    output:
        temp(outdir + "/samples/{sample}/{sample}.markdup.bam.bai")
    message:
        "Indexing duplicate-marked alignments: {wildcards.sample}"
    shell:
        "samtools index {input}"

rule bxstats_report:
    input:
        outdir + "/reports/BXstats/data/{sample}.bxstats.gz"
    output:	
        outdir + "/reports/BXstats/{sample}.bxstats.html"
    params:
        molecule_distance
    conda:
        os.getcwd() + "/.harpy_envs/r-env.yaml"
    message: 
        "Summarizing barcoded alignments: {wildcards.sample}"
    script:
        "report/BxStats.Rmd"

rule assign_molecules:
    input:
        bam = outdir + "/samples/{sample}/{sample}.markdup.bam",
        bai = outdir + "/samples/{sample}/{sample}.markdup.bam.bai"
    output:
        bam = outdir + "/{sample}.bam",
        bai = outdir + "/{sample}.bam.bai"
    params:
        molecule_distance
    conda:
        os.getcwd() + "/.harpy_envs/qc.yaml"
    message:
        "Assigning barcodes to molecules: {wildcards.sample}"
    script:
        "scripts/assignMI.py"

rule alignment_bxstats:
    input:
        bam = outdir + "/{sample}.bam",
        bai = outdir + "/{sample}.bam.bai"
    output: 
        outdir + "/reports/BXstats/data/{sample}.bxstats.gz"
    params:
        sample = lambda wc: d[wc.sample]
    conda:
        os.getcwd() + "/.harpy_envs/qc.yaml"
    message:
        "Calculating barcode alignment statistics: {wildcards.sample}"
    script:
        "scripts/bxStats.py"

rule alignment_coverage:
    input: 
        bed = f"Genome/{bn}.bed",
        bam = outdir + "/{sample}.bam"
    output: 
        outdir + "/reports/coverage/data/{sample}.cov.gz"
    threads: 
        2
    message:
        "Calculating genomic coverage: {wildcards.sample}"
    shell:
        "samtools bedcov -c {input} | gzip > {output}"

rule coverage_report:
    input:
        outdir + "/reports/coverage/data/{sample}.cov.gz"
    output:
        outdir + "/reports/coverage/{sample}.cov.html"
    conda:
        os.getcwd() + "/.harpy_envs/r-env.yaml"
    message:
        "Summarizing alignment coverage: {wildcards.sample}"
    script:
        "report/Gencov.Rmd"
    
rule general_alignment_stats:
    input:
        bam      = outdir + "/{sample}.bam",
        bai      = outdir + "/{sample}.bam.bai"
    output: 
        stats    = temp(outdir + "/reports/samtools_stats/{sample}.stats"),
        flagstat = temp(outdir + "/reports/samtools_flagstat/{sample}.flagstat")
    message:
        "Calculating alignment stats: {wildcards.sample}"
    shell:
        """
        samtools stats {input.bam} > {output.stats}
        samtools flagstat {input.bam} > {output.flagstat}
        """

rule samtools_reports:
    input: 
        expand(outdir + "/reports/samtools_{ext}/{sample}.{ext}", sample = samplenames, ext = ["stats", "flagstat"])
    output: 
        outdir + "/reports/minimap.stats.html"
    params:
        outdir
    conda:
        os.getcwd() + "/.harpy_envs/qc.yaml"
    message:
        "Summarizing samtools stats and flagstat"
    shell:
        """
        multiqc {params}/reports/samtools_stats {params}/reports/samtools_flagstat --force --quiet --title "Basic Alignment Statistics" --comment "This report aggregates samtools stats and samtools flagstats results for all alignments." --no-data-dir --filename {output} 2> /dev/null
        """

rule log_workflow:
    default_target: True
    input: 
        bams = expand(outdir + "/{sample}.{ext}", sample = samplenames, ext = ["bam","bam.bai"]),
        cov_reports = expand(outdir + "/reports/coverage/{sample}.cov.html", sample = samplenames) if not skipreports else [],
        bx_reports = expand(outdir + "/reports/BXstats/{sample}.bxstats.html", sample = samplenames) if not skipreports else []
    output:
        outdir + "/workflow/align.minimap.summary"
    params:
        quality = config["quality"],
        extra   = extra
    message:
        "Summarizing the workflow: {output}"
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy align module ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write(f"The directory with sequences: {seq_dir}\n\n")
            _ = f.write("Sequencing were aligned with Minimap2 using:\n")
            _ = f.write("    minimap2 -y " + " ".join([str(i) for i in params]) + " --sam-hit-only -R \"@RG\\tID:SAMPLE\\tSM:SAMPLE\" genome.mmi forward_reads reverse_reads |\n")
            _ = f.write("    samtools view -h -F 4 -q " + str(config["quality"]) + " |\n")
            _ = f.write("    samtools sort -T SAMPLE --reference genome -m 4G\n")
            _ = f.write("Duplicates in the alignments were marked following:\n")
            _ = f.write("    samtools collate \n")
            _ = f.write("    samtools fixmate\n")
            _ = f.write("    samtools sort \n")
            _ = f.write("    samtools markdup -S --barcode-tag BX\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")