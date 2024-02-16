import os
import re
import glob
import sys
from rich.panel import Panel
from rich import print as rprint

outdir      = "Align/bwa"
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
            title = "[bold]harpy align bwa",
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
            title = "[bold]harpy align bwa",
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

rule genome_bwa_index:
    input: 
        f"Genome/{bn}"
    output: 
        multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    log:
        f"Genome/{bn}.idx.log"
    conda:
        os.getcwd() + "/harpyenvs/align.yaml"
    message:
        "Indexing {input}"
    shell: 
        "bwa index {input} 2> {log}"

rule genome_make_windows:
    input:
        f"Genome/{bn}.fai"
    output: 
        f"Genome/{bn}.bed"
    message: 
        "Creating BED intervals from {input}"
    shell: 
        "makewindows.py -i {input} -w 10000 -o {output}"

rule align:
    input:
        forward_reads = get_fq1,
        reverse_reads = get_fq2,
        genome 		  = f"Genome/{bn}",
        genome_samidx = f"Genome/{bn_idx}",
        genome_idx 	  = multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output:  
        pipe(outdir + "/samples/{sample}/{sample}.raw.sam")
    log:
        outdir + "/logs/{sample}.bwa.log"
    params: 
        samps = lambda wc: d[wc.get("sample")],
        extra = extra
    benchmark:
        ".Benchmark/Mapping/bwa/align.{sample}.txt"
    threads:
        min(10, workflow.cores) - 2
    conda:
        os.getcwd() + "/harpyenvs/align.yaml"
    message:
        "Aligning sequences: {wildcards.sample}"
    shell:
        """
        bwa mem -C -t {threads} {params.extra} -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" {input.genome} {input.forward_reads} {input.reverse_reads} > {output} 2> {log}
        """
 
rule quality_filter:
    input:
        outdir + "/samples/{sample}/{sample}.raw.sam"
    output:
        temp(outdir + "/samples/{sample}/{sample}.bwa.sam")
    params: 
        quality = config["quality"]
    message:
        "Quality filtering alignments: {wildcards.sample}"
    shell:
        "samtools view -h -F 4 -q {params.quality} {input} > {output}"

rule collate:
    input:
        outdir + "/samples/{sample}/{sample}.bwa.sam"
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
        genome_samidx = f"Genome/{bn_idx}",
        genome_idx 	  = multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output:
        bam = temp(outdir + "/samples/{sample}/{sample}.sort.bam"),
        bai = temp(outdir + "/samples/{sample}/{sample}.sort.bam.bai")
    log:
        outdir + "/logs/{sample}.bwa.sort.log"
    params: 
        quality = config["quality"],
        tmpdir = lambda wc: outdir + "/." + d[wc.sample]
    message:
        "Sorting and quality filtering alignments: {wildcards.sample}"
    shell:
        """
        samtools sort -T {params.tmpdir} --reference {input.genome} -O bam -l 0 -m 4G --write-index -o {output.bam}##idx##{output.bai} {input.sam} 2> {log}
        rm -rf {params.tmpdir}
        """

rule markduplicates:
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
        os.getcwd() + "/harpyenvs/r-env.yaml"
    message: 
        "Generating summary of barcode alignment: {wildcards.sample}"
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
    message:
        "Assigning barcodes to molecules: {wildcards.sample}"
    shell:
        "assignMI.py -c {params} -i {input.bam} -o {output.bam}"

rule alignment_bxstats:
    input:
        bam = outdir + "/{sample}.bam",
        bai = outdir + "/{sample}.bam.bai"
    output: 
        outdir + "/reports/BXstats/data/{sample}.bxstats.gz"
    params:
        sample = lambda wc: d[wc.sample]
    message:
        "Calculating barcode alignment statistics: {wildcards.sample}"
    shell:
        "bxStats.py {input.bam} | gzip > {output}"

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
        os.getcwd() + "/harpyenvs/r-env.yaml"
    message:
        "Summarizing alignment coverage: {wildcards.sample}"
    script:
        "report/BwaGencov.Rmd"
    
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
        outdir + "/reports/bwa.stats.html",
    message:
        "Summarizing samtools stats and flagstat"
    shell:
        """
        multiqc Align/bwa/reports/samtools_stats Align/bwa/reports/samtools_flagstat --force --quiet --title "Basic Alignment Statistics" --comment "This report aggregates samtools stats and samtools flagstats results for all alignments." --no-data-dir --filename {output} 2> /dev/null
        """

rule log_runtime:
    output:
        outdir + "/workflow/align.workflow.summary"
    params:
        quality = config["quality"],
        extra   = extra
    message:
        "Creating record of relevant runtime parameters: {output}"
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy align module ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write(f"The directory with sequences: {seq_dir}\n\n")
            _ = f.write("Sequencing were aligned with BWA using:\n")
            _ = f.write("    bwa mem -C " + " ".join([str(i) for i in params]) + " -R \"@RG\\tID:SAMPLE\\tSM:SAMPLE\" genome forward_reads reverse_reads |\n")
            _ = f.write("    samtools view -h -F 4 -q " + str(config["quality"]) + " |\n")
            _ = f.write("    samtools sort -T SAMPLE --reference genome -m 4G\n")
            _ = f.write("Duplicates in the alignments were marked using sambamba:\n")
            _ = f.write("    sambamba markdup -l 0\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")

# conditionally add the reports to the output
results = list()
results.append(expand(outdir + "/{sample}.bam", sample = samplenames))
results.append(expand(outdir + "/{sample}.bam.bai", sample = samplenames))
results.append(outdir + "/workflow/align.workflow.summary")

if not skipreports:
    results.append(expand(outdir + "/reports/coverage/{sample}.cov.html", sample = samplenames))
    results.append(expand(outdir + "/reports/BXstats/{sample}.bxstats.html", sample = samplenames))
    results.append(outdir + "/reports/bwa.stats.html")

rule all:
    default_target: True
    input: 
        results
    message:
        "Checking for expected workflow output"
