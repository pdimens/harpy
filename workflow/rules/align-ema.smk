import os
import re
import glob
from pathlib import Path
from rich.panel import Panel
from rich import print as rprint

outdir      = "Align/ema"
seq_dir 	= config["seq_directory"]
nbins 		= config["EMA_bins"]
binrange    = ["%03d" % i for i in range(nbins)]
genomefile 	= config["genomefile"]
samplenames = config["samplenames"]
platform    = config["platform"]
whitelist   = config.get("whitelist", "") 
extra 		= config.get("extra", "") 
bn 			= os.path.basename(genomefile)
genome_zip  = True if bn.lower().endswith(".gz") else False
bn_idx      = f"{bn}.gzi" if genome_zip else f"{bn}.fai"

d = dict(zip(samplenames, samplenames))

def get_fq1(wildcards):
    # code that returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    lst = glob.glob(seq_dir + "/" + wildcards.sample + "*")
    return lst[0]

def get_fq2(wildcards):
    # code that returns a list of fastq files for read 2 based on *wildcards.sample*, e.g.
    lst = sorted(glob.glob(seq_dir + "/" + wildcards.sample + "*"))
    return lst[1]

onerror:
    print("")
    rprint(
        Panel(
            f"The workflow has terminated due to an error. See the log file below for more details.",
            title = "[bold]harpy align ema",
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
            title = "[bold]harpy align ema",
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

rule interleave:
    input:
        fw_reads = get_fq1,
        rv_reads = get_fq2
    output:
        pipe(outdir + "/.interleave/{sample}.interleave.fq")
    message:
        "Interleaving input fastq files: {wildcards.sample}"
    shell:
        "seqtk mergepe {input} > {output}"

use rule interleave as interleave2 with:
    output:
        outdir + "/.interleave/{sample}.interleave2.fq"

rule beadtag_count:
    input:
        outdir + "/.interleave/{sample}.interleave.fq"
    output: 
        counts = temp(outdir + "/bxcount/{sample}.ema-ncnt"),
        logs   = temp(outdir + "/logs/count/{sample}.count")
    params:
        prefix = lambda wc: outdir + "/bxcount/" + wc.get("sample"),
        beadtech = "-p" if platform == "haplotag" else f"-w {whitelist}",
        logdir = f"{outdir}/logs/count/"
    message:
        "Counting barcode frequency: {wildcards.sample}"
    conda:
        os.getcwd() + "/harpyenvs/align.yaml"
    shell:
        """
        mkdir -p {params.prefix} {params.logdir}
        ema count {params.beadtech} -o {params.prefix} < {input} 2> {output.logs}
        """

rule beadtag_summary:
    input: 
        countlog = expand(outdir + "/logs/count/{sample}.count", sample = samplenames)
    output:
        outdir + "/reports/reads.bxcounts.html"
    conda:
        os.getcwd() + "/harpyenvs/r-env.yaml"
    message:
        "Creating sample barcode validation report"
    script:
        "report/EmaCount.Rmd"

rule preprocess:
    input: 
        reads = outdir + "/.interleave/{sample}.interleave2.fq",
        emacounts  = outdir + "/bxcount/{sample}.ema-ncnt"
    output: 
        bins       = temp(expand(outdir + "/preproc/{{sample}}/ema-bin-{bin}", bin = binrange)),
        unbarcoded = temp(outdir + "/preproc/{sample}/ema-nobc")
    log:
        outdir + "/logs/preproc/{sample}.preproc.log"
    params:
        outdir = lambda wc: outdir + "/preproc/" + wc.get("sample"),
        bxtype = "-p" if platform == "haplotag" else f"-w {whitelist}",
        bins   = nbins
    threads:
        2
    conda:
        os.getcwd() + "/harpyenvs/align.yaml"
    message:
        "Preprocessing for EMA mapping: {wildcards.sample}"
    shell:
        """
        ema preproc {params.bxtype} -n {params.bins} -t {threads} -o {params.outdir} {input.emacounts} < {input.reads} 2>&1 |
            cat - > {log}
        """

rule align:
    input:
        readbin    = expand(outdir + "/preproc/{{sample}}/ema-bin-{bin}", bin = binrange),
        genome 	   = f"Genome/{bn}",
        geno_faidx = f"Genome/{bn_idx}",
        geno_idx   = multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output:
        pipe(outdir + "/align/{sample}.bc.raw.sam"),
    log:
        outdir + "/logs/{sample}.ema.align.log",
    params: 
        bxtype = f"-p {platform}",
        extra = extra
    threads:
        min(10, workflow.cores) - 2
    conda:
        os.getcwd() + "/harpyenvs/align.yaml"
    message:
        "Aligning barcoded sequences: {wildcards.sample}"
    shell:
        """
        ema align -t {threads} {params.extra} -d {params.bxtype} -r {input.genome} -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" -x {input.readbin} 2> {log}
        """

rule sort_raw_ema:
    input:
        aln = outdir + "/align/{sample}.bc.raw.sam",
        genome 	   = f"Genome/{bn}",
        geno_faidx = f"Genome/{bn_idx}",
        geno_idx   = multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output:
        aln = temp(outdir + "/align/{sample}.bc.bam"),
        idx = temp(outdir + "/align/{sample}.bc.bam.bai")
    log:
        outdir + "/logs/{sample}.ema.sort.log"
    params: 
        quality = config["quality"],
        tmpdir = lambda wc: outdir + "/align/." + d[wc.sample],
        extra = extra
    threads:
        2
    message:
        "Sorting and quality filtering alignments: {wildcards.sample}"
    shell:
        """
        samtools view -h -F 4 -q {params.quality} - | 
            samtools sort -T {params.tmpdir} --reference {input.genome} -O bam --write-index -m 4G -o {output.aln}##idx##{output.idx} - 2> {log}
        rm -rf {params.tmpdir}
        """

rule align_nobarcode:
    input:
        reads      = outdir + "/preproc/{sample}/ema-nobc",
        genome 	   = f"Genome/{bn}",
        geno_faidx = f"Genome/{bn_idx}",
        geno_idx   = multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output: 
        pipe(outdir + "/align/{sample}.nobc.raw.sam")
    log:
        outdir + "/logs/{sample}.bwa.align.log"
    benchmark:
        ".Benchmark/Mapping/ema/bwaAlign.{sample}.txt"
    threads:
        min(10, workflow.cores) - 2
    conda:
        os.getcwd() + "/harpyenvs/align.yaml"
    message:
        "Aligning unbarcoded sequences: {wildcards.sample}"
    shell:
        """
        bwa mem -t {threads} -C -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" {input.genome} {input.reads} > {output} 2> {log}
        """
        
rule sort_raw_nobarcode:
    input:
        sam        = outdir + "/align/{sample}.nobc.raw.sam",
        genome 	   = f"Genome/{bn}",
        geno_faidx = f"Genome/{bn_idx}",
        geno_idx   = multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output: 
        temp(outdir + "/align/{sample}.nobc.bam")
    log:
        outdir + "/logs/{sample}.bwa.sort.log"
    benchmark:
        ".Benchmark/Mapping/ema/bwaAlign.{sample}.txt"
    params:
        quality = config["quality"]
    threads:
        2
    message:
        "Aligning unbarcoded sequences: {wildcards.sample}"
    shell:        
        """
        samtools view -h -F 4 -q {params.quality} {input.sam} | 
            samtools sort -O bam -m 4G --reference {input.genome} -o {output} 2> {log}
        """

rule mark_duplicates:
    input:
        bam      = outdir + "/align/{sample}.nobc.bam"
    output: 
        bam      = temp(outdir + "/align/markdup/{sample}.markdup.nobc.bam"),
        bai      = temp(outdir + "/align/markdup/{sample}.markdup.nobc.bam.bai")
    log: 
        mdlog    = outdir + "/logs/markduplicates/{sample}.markdup.nobc.log",
        stats    = outdir + "/reports/samtools_stats/{sample}.nobc.stats",
        flagstat = outdir + "/reports/samtools_flagstat/{sample}.nobc.flagstat"
    conda:
        os.getcwd() + "/harpyenvs/align.yaml"
    threads:
        2
    message:
        "Marking duplicates in unbarcoded alignments: {wildcards.sample}"
    shell:
        """
        sambamba markdup -t {threads} -l 4 {input} {output.bam} 2> {log.mdlog}
        samtools stats {output.bam} > {log.stats}
        samtools flagstat {output.bam} > {log.flagstat}
        """

rule coverage_stats:
    input: 
        bed     = f"Genome/{bn}.bed",
        nobx    = outdir + "/align/markdup/{sample}.markdup.nobc.bam",
        nobxbai = outdir + "/align/markdup/{sample}.markdup.nobc.bam.bai",
        bx      = outdir + "/align/{sample}.bc.bam",
        bxbai   = outdir + "/align/{sample}.bc.bam.bai"
    output: 
        outdir + "/reports/coverage/data/{sample}.cov.gz"
    threads:
        2
    message:
        "Calculating genomic coverage: {wildcards.sample}"
    shell:
        "samtools bedcov -c {input.bed} {input.bx} {input.nobx} | gzip > {output}"

rule coverage_report:
    input: 
        outdir + "/reports/coverage/data/{sample}.cov.gz",
    output:
        outdir + "/reports/coverage/{sample}.cov.html"
    conda:
        os.getcwd() + "/harpyenvs/r-env.yaml"
    message:
        "Creating report of alignment coverage: {wildcards.sample}"
    script:
        "report/EmaGencov.Rmd"

rule merge_alignments:
    input:
        aln_bc   = outdir + "/align/{sample}.bc.bam",
        idx_bc   = outdir + "/align/{sample}.bc.bam.bai",
        aln_nobc = outdir + "/align/markdup/{sample}.markdup.nobc.bam",
        idx_nobc = outdir + "/align/markdup/{sample}.markdup.nobc.bam.bai"
    output: 
        bam 	 = temp(outdir + "/align/{sample}.unsort.bam"),
        bai 	 = temp(outdir + "/align/{sample}.unsort.bam.bai")
    conda:
        os.getcwd() + "/harpyenvs/align.yaml"
    threads:
        min(10, workflow.cores)
    message:
        "Merging all alignments: {wildcards.sample}"
    shell:
        "sambamba merge -t {threads} {output.bam} {input.aln_bc} {input.aln_nobc} 2> /dev/null"

rule sort_merge:
    input:
        bam    = outdir + "/align/{sample}.unsort.bam",
        genome = f"Genome/{bn}"
    output:
        bam = outdir + "/{sample}.bam",
        bai = outdir + "/{sample}.bam.bai"
    threads:
        2
    message:
        "Sorting merged barcoded alignments: {wildcards.sample}"
    shell:
        "samtools sort -@ {threads} -O bam --reference {input.genome} -m 4G --write-index -o {output.bam}##idx##{output.bai} {input.bam} 2> /dev/null"

rule bx_stats_alignments:
    input:
        bam = outdir + "/{sample}.bam",
        bai = outdir + "/{sample}.bam.bai"
    output: 
        outdir + "/reports/BXstats/data/{sample}.bxstats.gz"
    message:
        "Calculating barcode alignment statistics: {wildcards.sample}"
    shell:
        "bxStats.py {input.bam} | gzip > {output}"

rule bx_stats_report:
    input:
        outdir + "/reports/BXstats/data/{sample}.bxstats.gz"
    output:	
        outdir + "/reports/BXstats/{sample}.bxstats.html"
    params:
        "none"
    conda:
        os.getcwd() + "/harpyenvs/r-env.yaml"
    message: 
        "Generating summary of barcode alignment: {wildcards.sample}"
    script:
        "report/BxStats.Rmd"

rule general_stats:
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

rule collate_samtools_stats:
    input: 
        expand(outdir + "/reports/samtools_{ext}/{sample}.{ext}", sample = samplenames, ext = ["stats", "flagstat"]),
    output: 
        outdir + "/reports/ema.stats.html"
    message:
        "Summarizing samtools stats and flagstat"
    shell:
        """
        multiqc Align/ema/reports/samtools_stats Align/ema/reports/samtools_flagstat --force --quiet --title "Basic Alignment Statistics" --comment "This report aggregates samtools stats and samtools flagstats results for all alignments." --no-data-dir --filename {output} 2> /dev/null
        """

rule log_runtime:
    output:
        outdir + "/workflow/align.workflow.summary"
    params:
        beadtech = "-p" if platform == "haplotag" else f"-w {whitelist}"
    message:
        "Creating record of relevant runtime parameters: {output}"
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy align module ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write(f"The directory with sequences: {seq_dir}\n")
            _ = f.write("Barcodes were counted and validated with EMA using:\n")
            _ = f.write(f"    seqtk mergepe forward.fq.gz reverse.fq.gz | ema count {params.beadtech}\n")
            _ = f.write("Barcoded sequences were binned with EMA using:\n")
            _ = f.write(f"    seqtk mergepe forward.fq.gz reverse.fq.gz | ema preproc {params.beadtech} -n {nbins}\n")
            _ = f.write("Barcoded bins were aligned with ema align using:\n")
            _ = f.write(f"    ema align " + extra + " -d -p " + platform + " -R \"@RG\\tID:SAMPLE\\tSM:SAMPLE\" |\n")
            _ = f.write("    samtools view -h -F 4 -q " + str(config["quality"]) + " - |\n") 
            _ = f.write("    samtools sort --reference genome -m 4G\n\n")
            _ = f.write("Invalid/non barcoded sequences were aligned with BWA using:\n")
            _ = f.write("    bwa mem -C -R \"@RG\\tID:SAMPLE\\tSM:SAMPLE\" genome forward_reads reverse_reads\n")
            _ = f.write("Duplicats in non-barcoded alignments were marked using:\n")
            _ = f.write("    sambamba markdup -l 0\n")
            _ = f.write("Alignments were merged using:\n")
            _ = f.write("    sambamba merge output.bam barcode.bam nobarcode.bam\n")
            _ = f.write("Merged alignments were sorted using:\n")
            _ = f.write("    samtools sort merged.bam\n")

rule all:
    default_target: True
    input: 
        bam = expand(outdir + "/{sample}.bam", sample = samplenames),
        bai = expand(outdir + "/{sample}.bam.bai", sample = samplenames),
        covstats = expand(outdir + "/reports/coverage/{sample}.cov.html", sample = samplenames),
        bxstats = expand(outdir + "/reports/BXstats/{sample}.bxstats.html", sample = samplenames),
        bxcounts = f"{outdir}/reports/reads.bxcounts.html",
        emastats = f"{outdir}/reports/ema.stats.html",
        runlog = f"{outdir}/workflow/align.workflow.summary"
    message:
        "Checking for expected workflow output"