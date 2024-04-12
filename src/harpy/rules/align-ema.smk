import os
import re
import glob
from pathlib import Path
from rich.panel import Panel
from rich import print as rprint

outdir      = config["output_directory"]
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
skipreports = config["skipreports"]

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
        os.getcwd() + "/.harpy_envs/align.yaml"
    message:
        "Indexing {input}"
    shell: 
        "bwa index {input} 2> {log}"

rule genome_make_windows:
    input:
        f"Genome/{bn}"
    output: 
        f"Genome/{bn}.bed"
    message: 
        "Creating BED intervals from {input}"
    shell: 
        "makeWindows.py -i {input} -w 10000 -o {output}"

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
        os.getcwd() + "/.harpy_envs/align.yaml"
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
        os.getcwd() + "/.harpy_envs/r-env.yaml"
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
        os.getcwd() + "/.harpy_envs/align.yaml"
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
        os.getcwd() + "/.harpy_envs/align.yaml"
    message:
        "Aligning barcoded sequences: {wildcards.sample}"
    shell:
        """
        ema align -t {threads} {params.extra} -d {params.bxtype} -r {input.genome} -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" -x {input.readbin} > {output} 2> {log}
        """

rule sort_raw_ema:
    input:
        sam        = outdir + "/align/{sample}.bc.raw.sam",
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
        samtools view -h -F 4 -q {params.quality} {input.sam} | 
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
        os.getcwd() + "/.harpy_envs/align.yaml"
    message:
        "Aligning unbarcoded sequences: {wildcards.sample}"
    shell:
        """
        bwa mem -t {threads} -C -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" {input.genome} {input.reads} > {output} 2> {log}
        """

rule quality_filter:
    input:
        outdir + "/align/{sample}.nobc.raw.sam"
    output:
        temp(outdir + "/align/{sample}.bwa.nobc.sam")
    params: 
        quality = config["quality"]
    message:
        "Quality filtering alignments: {wildcards.sample}"
    shell:
        "samtools view -h -F 4 -q {params.quality} {input} > {output}"

rule collate:
    input:
        outdir + "/align/{sample}.bwa.nobc.sam"
    output:
        temp(outdir + "/align/{sample}.collate.nobc.bam")
    message:
        "Collating alignments: {wildcards.sample}"
    shell:
        "samtools collate -o {output} {input} 2> /dev/null"

rule fix_mates:
    input:
        outdir + "/align/{sample}.collate.nobc.bam"
    output:
        temp(outdir + "/align/{sample}.fixmate.nobc.bam")
    message:
        "Fixing mates in alignments: {wildcards.sample}"
    shell:
        "samtools fixmate -m {input} {output} 2> /dev/null"

rule sort_nobc_alignments:
    input:
        sam           = outdir + "/align/{sample}.fixmate.nobc.bam",
        genome 		  = f"Genome/{bn}",
        genome_samidx = f"Genome/{bn_idx}",
        genome_idx 	  = multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output:
        bam = temp(outdir + "/align/{sample}.sort.nobc.bam"),
        bai = temp(outdir + "/align/{sample}.sort.nobc.bam.bai")
    log:
        outdir + "/logs/{sample}.bwa.sort.log"
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
        bam = outdir + "/align/{sample}.sort.nobc.bam",
        bai = outdir + "/align/{sample}.sort.nobc.bam.bai"
    output:
        temp(outdir + "/align/{sample}.markdup.nobc.bam")
    log:
        outdir + "/logs/{sample}.markdup.log"
    threads:
        2
    message:
        "Marking duplicates in alignments alignment: {wildcards.sample}"
    shell:
        "samtools markdup -@ {threads} -S -f {log} {input.bam} {output}  2> /dev/null"

rule index_markdups:
    input:
        outdir + "/align/{sample}.markdup.nobc.bam"
    output:
        temp(outdir + "/align/{sample}.markdup.nobc.bam.bai")
    message:
        "Indexing duplicate-marked alignments: {wildcards.sample}"
    shell:
        "samtools index {input}"

rule coverage_stats:
    input: 
        bed     = f"Genome/{bn}.bed",
        nobx    = outdir + "/align/{sample}.markdup.nobc.bam",
        nobxbai = outdir + "/align/{sample}.markdup.nobc.bam.bai",
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
        os.getcwd() + "/.harpy_envs/r-env.yaml"
    message:
        "Creating report of alignment coverage: {wildcards.sample}"
    script:
        "report/EmaGencov.Rmd"

rule concatenate_alignments:
    input:
        aln_bc   = outdir + "/align/{sample}.bc.bam",
        idx_bc   = outdir + "/align/{sample}.bc.bam.bai",
        aln_nobc = outdir + "/align/{sample}.markdup.nobc.bam",
        idx_nobc = outdir + "/align/{sample}.markdup.nobc.bam.bai"
    output: 
        bam 	 = temp(outdir + "/align/{sample}.concat.unsort.bam")
    threads:
        2
    message:
        "Concatenating barcoded and unbarcoded alignments: {wildcards.sample}"
    shell:
        "samtools cat -@ {threads} {input.aln_bc} {input.aln_nobc} > {output.bam}"

rule sort_concatenated:
    input:
        bam    = outdir + "/align/{sample}.concat.unsort.bam",
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
    conda:
        os.getcwd() + "/.harpy_envs/qc.yaml"
    message:
        "Calculating barcode alignment statistics: {wildcards.sample}"
    script:
        "scripts/bxStats.py"

rule bx_stats_report:
    input:
        outdir + "/reports/BXstats/data/{sample}.bxstats.gz"
    output:	
        outdir + "/reports/BXstats/{sample}.bxstats.html"
    params:
        "none"
    conda:
        os.getcwd() + "/.harpy_envs/r-env.yaml"
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
    params:
        outdir
    conda:
        os.getcwd() + "/.harpy_envs/qc.yaml"
    message:
        "Summarizing samtools stats and flagstat"
    shell:
        """
        multiqc {outdir}/reports/samtools_stats {outdir}/reports/samtools_flagstat --force --quiet --title "Basic Alignment Statistics" --comment "This report aggregates samtools stats and samtools flagstats results for all alignments." --no-data-dir --filename {output} 2> /dev/null
        """

# conditionally add the reports to the output

rule log_workflow:
    default_target: True
    input:
        bams = expand(outdir + "/{sample}.{ext}", sample = samplenames, ext = [ "bam", "bam.bai"] ),
        cov_report = expand(outdir + "/reports/coverage/{sample}.cov.html", sample = samplenames) if not skipreports else [],
        bx_report = expand(outdir + "/reports/BXstats/{sample}.bxstats.html", sample = samplenames) if not skipreports else [],
        bx_counts = f"{outdir}/reports/reads.bxcounts.html" if not skipreports else [],
        agg_report = f"{outdir}/reports/ema.stats.html" if not skipreports else []
    output:
        outdir + "/workflow/align.ema.summary"
    params:
        beadtech = "-p" if platform == "haplotag" else f"-w {whitelist}"
    message:
        "Summarizing the workflow: {output}"
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
            _ = f.write("Duplicates in non-barcoded alignments were marked following:\n")
            _ = f.write("    samtools collate \n")
            _ = f.write("    samtools fixmate\n")
            _ = f.write("    samtools sort \n")
            _ = f.write("    samtools markdup -S\n")
            _ = f.write("Alignments were merged using:\n")
            _ = f.write("    samtools cat barcode.bam nobarcode.bam > concat.bam\n")
            _ = f.write("Merged alignments were sorted using:\n")
            _ = f.write("    samtools sort concat.bam\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")