containerized: "docker://pdimens/harpy:latest"

import os
import re
import logging

onstart:
    logger.logger.addHandler(logging.FileHandler(config["snakemake_log"]))
onsuccess:
    os.remove(logger.logfile)
onerror:
    os.remove(logger.logfile)
wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

outdir      = config["output_directory"]
envdir      = os.getcwd() + "/.harpy_envs"
genomefile 	= config["inputs"]["genome"]
fqlist      = config["inputs"]["fastq"]
extra 		= config.get("extra", "") 
bn 			= os.path.basename(genomefile)
if bn.lower().endswith(".gz"):
    bn = bn[:-3]
skip_reports = config["skip_reports"]
windowsize  = config["depth_windowsize"]
molecule_distance = config["molecule_distance"]
keep_unmapped = config["keep_unmapped"]
readlen = config["average_read_length"]
autolen = isinstance(readlen, str)
bn_r = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
samplenames = {re.sub(bn_r, "", os.path.basename(i), flags = re.IGNORECASE) for i in fqlist}
d = dict(zip(samplenames, samplenames))

def get_fq(wildcards):
    # returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    r = re.compile(fr".*/({re.escape(wildcards.sample)}){bn_r}", flags = re.IGNORECASE)
    return sorted(list(filter(r.match, fqlist))[:2])

rule process_genome:
    input:
        genomefile
    output: 
        f"Genome/{bn}"
    container:
        None
    shell: 
        "seqtk seq {input} > {output}"

rule index_genome:
    input: 
        f"Genome/{bn}"
    output: 
        f"Genome/{bn}.fai",
    log:
        f"Genome/{bn}.faidx.log"
    container:
        None
    shell: 
        "samtools faidx --fai-idx {output} {input} 2> {log}"

rule strobe_index:
    input: 
        f"Genome/{bn}"
    output:
        f"Genome/{bn}.r{readlen}.sti"
    log:
        f"Genome/{bn}.r{readlen}.sti.log"
    params:
        readlen
    conda:
        f"{envdir}/align.yaml"
    threads:
        2
    shell: 
        "strobealign --create-index -t {threads} -r {params} {input} 2> {log}"

rule align:
    input:
        fastq = get_fq,
        genome   = f"Genome/{bn}",
        genome_index   = f"Genome/{bn}.r{readlen}.sti" if not autolen else []
    output:  
        temp(outdir + "/samples/{sample}/{sample}.strobe.sam")
    log:
        outdir + "/logs/strobealign/{sample}.strobealign.log"
    params: 
        samps = lambda wc: d[wc.get("sample")],
        readlen = "" if autolen else f"--use-index -r {readlen}",
        quality = config["quality"],
        unmapped_strobe = "" if keep_unmapped else "-U",
        unmapped = "" if keep_unmapped else "-F 4",
        extra = extra
    benchmark:
        ".Benchmark/Mapping/strobealign/align.{sample}.txt"
    threads:
        max(5, workflow.cores - 1)
    conda:
        f"{envdir}/align.yaml"
    shell:
        """
        strobealign {params.readlen} -N 2 -t {threads} {params.unmapped_strobe} -C --rg-id={wildcards.sample} --rg=SM:{wildcards.sample} {params.extra} {input.genome} {input.fastq} 2> {log} |
            samtools view -h {params.unmapped} -q {params.quality} > {output} 
        """

rule mark_duplicates:
    input:
        sam    = outdir + "/samples/{sample}/{sample}.strobe.sam",
        genome = f"Genome/{bn}",
        faidx  = f"Genome/{bn}.fai"
    output:
        temp(outdir + "/samples/{sample}/{sample}.markdup.bam")
    log:
        outdir + "/logs/markdup/{sample}.markdup.log"
    params: 
        tmpdir = lambda wc: outdir + "/." + d[wc.sample]
    resources:
        mem_mb = 2000
    threads:
        2
    container:
        None
    shell:
        """
        samtools collate -O -u {input.sam} |
            samtools fixmate -m -u - - |
            samtools sort -T {params.tmpdir} -u --reference {input.genome} -l 0 -m {resources.mem_mb}M - |
            samtools markdup -@ {threads} -S --barcode-tag BX -f {log} - {output}
        rm -rf {params.tmpdir}
        """

rule index_duplicates:
    input:
        outdir + "/samples/{sample}/{sample}.markdup.bam"
    output:
        temp(outdir + "/samples/{sample}/{sample}.markdup.bam.bai")
    container:
        None
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
    shell:
        "assign_mi.py -o {output.bam} -c {params} {input.bam}"

rule barcode_stats:
    input:
        bam = outdir + "/{sample}.bam",
        bai = outdir + "/{sample}.bam.bai"
    output: 
        outdir + "/reports/data/bxstats/{sample}.bxstats.gz"
    params:
        sample = lambda wc: d[wc.sample]
    container:
        None
    shell:
        "bx_stats.py -o {output} {input.bam}"

rule calculate_depth:
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

rule sample_reports:
    input:
        outdir + "/reports/data/bxstats/{sample}.bxstats.gz",
        outdir + "/reports/data/coverage/{sample}.cov.gz"
    output:	
        outdir + "/reports/{sample}.html"
    log:
        logfile = outdir + "/logs/reports/{sample}.alignstats.log"
    params:
        molecule_distance
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/align_stats.Rmd"

rule general_stats:
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

rule samtools_report:
    input: 
        collect(outdir + "/reports/data/samtools_{ext}/{sample}.{ext}", sample = samplenames, ext = ["stats", "flagstat"])
    output: 
        outdir + "/reports/strobealign.stats.html"
    params:
        outdir = f"{outdir}/reports/data/samtools_stats {outdir}/reports/data/samtools_flagstat",
        options = "--no-version-check --force --quiet --no-data-dir",
        title = "--title \"Basic Alignment Statistics\"",
        comment = "--comment \"This report aggregates samtools stats and samtools flagstats results for all alignments. Samtools stats ignores alignments marked as duplicates.\""
    conda:
        f"{envdir}/qc.yaml"
    shell:
        "multiqc  {params} --filename {output} 2> /dev/null"

rule barcode_report:
    input:
        collect(outdir + "/reports/data/bxstats/{sample}.bxstats.gz", sample = samplenames)
    output:	
        outdir + "/reports/barcodes.summary.html"
    log:
        logfile = outdir + "/logs/reports/bxstats.report.log"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/align_bxstats.Rmd"

rule workflow_summary:
    default_target: True
    input: 
        bams = collect(outdir + "/{sample}.{ext}", sample = samplenames, ext = ["bam","bam.bai"]),
        samtools =  outdir + "/reports/strobealign.stats.html" if not skip_reports else [] ,
        reports = collect(outdir + "/reports/{sample}.html", sample = samplenames) if not skip_reports else [],
        bx_report = outdir + "/reports/barcodes.summary.html" if (not skip_reports or len(samplenames) == 1) else []
    params:
        readlen = readlen,
        quality = config["quality"],
        unmapped_strobe = "" if keep_unmapped else "-U",
        unmapped = "" if keep_unmapped else "-F 4",
        extra   = extra
    run:
        if autolen:
            strobetext = f"Sequences were aligned with strobealign using:\n"
            strobetext += f"    strobealign -U -C --rg-id=SAMPLE --rg=SM:SAMPLE {params.extra} genome reads.F.fq reads.R.fq |\n"
        else:
            strobetext = "The genome index was created using:\n"
            strobetext += f"    strobealign --create-index -r {params.readlen} genome\n\n"
            strobetext += "Sequences were aligned with strobealign using:\n"
            strobetext += f"    strobealign --use-index {params.unmapped_strobe} -N 2 -C --rg=SM:SAMPLE {params.extra} genome reads.F.fq reads.R.fq |\n"
        summary_template = f"""
The harpy align strobealign workflow ran using these parameters:

The provided genome: {bn}

{strobetext}
    samtools view -h {params.unmapped} -q {params.quality}

Duplicates in the alignments were marked following:
    samtools collate |
    samtools fixmate |
    samtools sort -m 2000M |
    samtools markdup -S --barcode-tag BX

The Snakemake workflow was called via command line:
    {config["workflow_call"]}
"""
        with open(outdir + "/workflow/align.strobealign.summary", "w") as f:
            f.write(summary_template)