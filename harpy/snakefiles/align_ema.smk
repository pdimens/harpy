containerized: "docker://pdimens/harpy:latest"

import os
import re
import shutil
import logging as pylogging
from pathlib import Path

outdir      = config["output_directory"]
nbins 		= config["EMA_bins"]
binrange    = ["%03d" % i for i in range(nbins)]
fqlist       = config["inputs"]["fastq"]
genomefile 	= config["inputs"]["genome"]
platform    = config["platform"]
whitelist   = config["inputs"].get("whitelist", "") 
extra 		= config.get("extra", "") 
bn 			= os.path.basename(genomefile)
genome_zip  = True if bn.lower().endswith(".gz") else False
bn_idx      = f"{bn}.gzi" if genome_zip else f"{bn}.fai"
envdir      = os.getcwd() + "/.harpy_envs"
windowsize  = config["depth_windowsize"]
keep_unmapped = config["keep_unmapped"]
skipreports = config["skip_reports"]
snakemake_log = config["snakemake_log"]

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

onstart:
    extra_logfile_handler = pylogging.FileHandler(snakemake_log)
    logger.logger.addHandler(extra_logfile_handler)

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
    shell: 
        """
        if (file {input} | grep -q compressed ) ;then
            # is regular gzipped, needs to be BGzipped
            seqtk seq {input} | bgzip -c > {output}
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
    shell: 
        """
        if [ "{params}" = "True" ]; then
            samtools faidx --gzi-idx {output.gzi} --fai-idx {output.fai} {input} 2> {log}
        else
            samtools faidx --fai-idx {output.fai} {input} 2> {log}
        fi
        """

rule genome_bwa_index:
    input: 
        f"Genome/{bn}"
    output: 
        multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    log:
        f"Genome/{bn}.idx.log"
    conda:
        f"{envdir}/align.yaml"
    shell: 
        "bwa index {input} 2> {log}"

rule ema_count:
    input:
        get_fq
    output: 
        counts = temp(outdir + "/ema_count/{sample}.ema-ncnt"),
        logs   = temp(outdir + "/logs/count/{sample}.count")
    params:
        prefix = lambda wc: outdir + "/ema_count/" + wc.get("sample"),
        beadtech = "-p" if platform == "haplotag" else f"-w {whitelist}",
        logdir = f"{outdir}/logs/ema_count/"
    conda:
        f"{envdir}/align.yaml"
    shell:
        """
        mkdir -p {params.prefix} {params.logdir}
        seqtk mergepe {input} |
            ema count {params.beadtech} -o {params.prefix} 2> {output.logs}
        """

rule ema_preprocess:
    input: 
        reads = get_fq,
        emacounts  = outdir + "/ema_count/{sample}.ema-ncnt"
    output: 
        bins       = temp(collect(outdir + "/ema_preproc/{{sample}}/ema-bin-{bin}", bin = binrange)),
        unbarcoded = temp(outdir + "/ema_preproc/{sample}/ema-nobc")
    log:
        outdir + "/logs/ema_preproc/{sample}.preproc.log"
    params:
        outdir = lambda wc: outdir + "/ema_preproc/" + wc.get("sample"),
        bxtype = "-p" if platform == "haplotag" else f"-w {whitelist}",
        bins   = nbins
    threads:
        2
    conda:
        f"{envdir}/align.yaml"
    shell:
        """
        seqtk mergepe {input.reads} |
            ema preproc {params.bxtype} -n {params.bins} -t {threads} -o {params.outdir} {input.emacounts} 2>&1 |
            cat - > {log}
        """

rule ema_align:
    input:
        readbin    = collect(outdir + "/ema_preproc/{{sample}}/ema-bin-{bin}", bin = binrange),
        genome 	   = f"Genome/{bn}",
        geno_faidx = f"Genome/{bn_idx}",
        geno_idx   = multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output:
        aln = temp(outdir + "/ema_align/{sample}.bc.bam"),
        idx = temp(outdir + "/ema_align/{sample}.bc.bam.bai")
    log:
        ema  = outdir + "/logs/align/{sample}.ema.align.log",
        sort = outdir + "/logs/align/{sample}.ema.sort.log",
    resources:
        mem_mb = 500
    params: 
        bxtype = f"-p {platform}",
        tmpdir = lambda wc: outdir + "/." + d[wc.sample],
        quality = config["quality"],
        unmapped = "" if keep_unmapped else "-F 4",
        extra = extra
    threads:
        10
    conda:
        f"{envdir}/align.yaml"
    shell:
        """
        ema align -t {threads} {params.extra} -d {params.bxtype} -r {input.genome} -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" -x {input.readbin} 2> {log.ema} |
            samtools view -h {params.unmapped} -q {params.quality} | 
            samtools sort -T {params.tmpdir} --reference {input.genome} -O bam --write-index -m {resources.mem_mb}M -o {output.aln}##idx##{output.idx} - 2> {log.sort}
        rm -rf {params.tmpdir}
        """

rule bwa_align:
    input:
        reads      = outdir + "/ema_preproc/{sample}/ema-nobc",
        genome 	   = f"Genome/{bn}",
        geno_faidx = f"Genome/{bn_idx}",
        geno_idx   = multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output: 
        temp(outdir + "/bwa_align/{sample}.bwa.nobc.sam")
    log:
        outdir + "/logs/align/{sample}.bwa.align.log"
    params:
        quality = config["quality"],
        unmapped = "" if keep_unmapped else "-F 4"
    benchmark:
        ".Benchmark/Mapping/ema/bwaAlign.{sample}.txt"
    threads:
        10
    conda:
        f"{envdir}/align.yaml"
    shell:
        """
        bwa mem -t {threads} -v2 -C -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" {input.genome} {input.reads} 2> {log} |
            samtools view -h {params.unmapped} -q {params.quality} > {output}
        """

rule bwa_markdups:
    input:
        sam    = outdir + "/bwa_align/{sample}.bwa.nobc.sam",
        genome = f"Genome/{bn}",
        faidx  = f"Genome/{bn_idx}"
    output:
        temp(outdir + "/bwa_align/{sample}.markdup.nobc.bam")
    log:
        outdir + "/logs/align/{sample}.markdup.log"
    params: 
        tmpdir = lambda wc: outdir + "/." + d[wc.sample]
    resources:
        mem_mb = 500
    container:
        None
    threads:
        2
    shell:
        """
        samtools collate -O -u {input.sam} |
            samtools fixmate -m -u - - |
            samtools sort -T {params.tmpdir} -u --reference {input.genome} -l 0 -m {resources.mem_mb}M - |
            samtools markdup -@ {threads} -S --barcode-tag BX -f {log} - {output}
        rm -rf {params.tmpdir}
        """

rule bwa_markdups_index:
    input:
        outdir + "/bwa_align/{sample}.markdup.nobc.bam"
    output:
        temp(outdir + "/bwa_align/{sample}.markdup.nobc.bam.bai")
    container:
        None
    shell:
        "samtools index {input}"

rule concatenate_alignments:
    input:
        aln_bc   = outdir + "/ema_align/{sample}.bc.bam",
        idx_bc   = outdir + "/ema_align/{sample}.bc.bam.bai",
        aln_nobc = outdir + "/bwa_align/{sample}.markdup.nobc.bam",
        idx_nobc = outdir + "/bwa_align/{sample}.markdup.nobc.bam.bai",
        genome   = f"Genome/{bn}"
    output: 
        bam = outdir + "/{sample}.bam",
        bai = outdir + "/{sample}.bam.bai"
    threads:
        2
    resources:
        mem_mb = 500
    container:
        None
    shell:
        """
        samtools cat -@ 1 {input.aln_bc} {input.aln_nobc} |
            samtools sort -@ 1 -O bam --reference {input.genome} -m {resources.mem_mb}M --write-index -o {output.bam}##idx##{output.bai} -
        """

rule coverage:
    input: 
        bam = outdir + "/{sample}.bam",
        bai = outdir + "/{sample}.bam.bai"
    output: 
        outdir + "/reports/data/coverage/{sample}.cov.gz"
    params:
        windowsize
    container:
        None
    shell:
        "samtools depth -a {input.bam} | depth_windows.py {params} | gzip > {output}"

rule bx_stats:
    input:
        bam = outdir + "/{sample}.bam",
        bai = outdir + "/{sample}.bam.bai"
    output: 
        outdir + "/reports/data/bxstats/{sample}.bxstats.gz"
    container:
        None
    shell:
        "bx_stats.py -o {output} {input.bam}"

rule report_persample:
    input:
        outdir + "/reports/data/bxstats/{sample}.bxstats.gz",
        outdir + "/reports/data/coverage/{sample}.cov.gz"
    output:	
        outdir + "/reports/{sample}.html"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/align_stats.Rmd"

rule stats:
    input: 		
        bam      = outdir + "/{sample}.bam",
        bai      = outdir + "/{sample}.bam.bai"
    output:
        stats    = temp(outdir + "/reports/data/samtools_stats/{sample}.stats"),
        flagstat = temp(outdir + "/reports/data/samtools_flagstat/{sample}.flagstat")
    container:
        None
    shell:
        """
        samtools stats -d {input.bam} > {output.stats}
        samtools flagstat {input.bam} > {output.flagstat}
        """

rule report_samtools:
    input: 
        collect(outdir + "/reports/data/samtools_{ext}/{sample}.{ext}", sample = samplenames, ext = ["stats", "flagstat"]),
    output: 
        outdir + "/reports/ema.stats.html"
    params:
        outdir = f"{outdir}/reports/data/samtools_stats {outdir}/reports/data/samtools_flagstat",
        options = "--no-version-check --force --quiet --no-data-dir",
        title = "--title \"Basic Alignment Statistics\"",
        comment = "--comment \"This report aggregates samtools stats and samtools flagstats results for all alignments. Samtools stats ignores alignments marked as duplicates.\""
    conda:
        f"{envdir}/qc.yaml"
    shell:
        "multiqc {params} --filename {output} 2> /dev/null"

rule report_bx:
    input:
        collect(outdir + "/reports/data/bxstats/{sample}.bxstats.gz", sample = samplenames)
    output:	
        outdir + "/reports/barcodes.summary.html"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/align_bxstats.Rmd"

rule workflow_summary:
    default_target: True
    input:
        bams = collect(outdir + "/{sample}.{ext}", sample = samplenames, ext = [ "bam", "bam.bai"] ),
        cov_report = collect(outdir + "/reports/{sample}.html", sample = samplenames) if not skipreports else [],
        agg_report = f"{outdir}/reports/ema.stats.html" if not skipreports else [],
        bx_report = outdir + "/reports/barcodes.summary.html" if (not skipreports or len(samplenames) == 1) else []
    params:
        beadtech = "-p" if platform == "haplotag" else f"-w {whitelist}",
        unmapped = "" if keep_unmapped else "-F 4"
    message:
        "Summarizing the workflow: {output}"
    run:
        with open(outdir + "/workflow/align.ema.summary", "w") as f:
            _ = f.write("The harpy align ema workflow ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write("Barcodes were counted and validated with EMA using:\n")
            _ = f.write(f"    seqtk mergepe forward.fq.gz reverse.fq.gz | ema count {params.beadtech}\n")
            _ = f.write("Barcoded sequences were binned with EMA using:\n")
            _ = f.write(f"    seqtk mergepe forward.fq.gz reverse.fq.gz | ema preproc {params.beadtech} -n {nbins}\n")
            _ = f.write("Barcoded bins were aligned with ema align using:\n")
            _ = f.write(f"    ema align " + extra + " -d -p " + platform + " -R \"@RG\\tID:SAMPLE\\tSM:SAMPLE\" |\n")
            _ = f.write("      samtools view -h {params.unmapped} -q " + str(config["quality"]) + " - |\n") 
            _ = f.write("      samtools sort --reference genome -m 2000M\n\n")
            _ = f.write("Invalid/non barcoded sequences were aligned with BWA using:\n")
            _ = f.write("    bwa mem -C -v2 -R \"@RG\\tID:SAMPLE\\tSM:SAMPLE\" genome forward_reads reverse_reads\n")
            _ = f.write("      samtools view -h {params.unmapped} -q " + str(config["quality"]) + " - |\n")
            _ = f.write("      samtools sort --reference genome -m 2000M\n\n") 
            _ = f.write("Duplicates in non-barcoded alignments were marked following:\n")
            _ = f.write("    samtools collate |\n")
            _ = f.write("      samtools fixmate |\n")
            _ = f.write("      samtools sort -m 2000M |\n")
            _ = f.write("      samtools markdup -S\n")
            _ = f.write("Alignments were merged using:\n")
            _ = f.write("    samtools cat barcode.bam nobarcode.bam > concat.bam\n")
            _ = f.write("Merged alignments were sorted using:\n")
            _ = f.write("    samtools sort -m 2000M concat.bam\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")