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
envdir      = os.path.join(os.getcwd(), ".harpy_envs")
genomefile 	= config["inputs"]["genome"]
fqlist       = config["inputs"]["fastq"]
molecule_distance = config["molecule_distance"]
keep_unmapped = config["keep_unmapped"]
extra 		= config.get("extra", "") 
bn 			= os.path.basename(genomefile)
genome_zip  = True if bn.lower().endswith(".gz") else False
bn_idx      = f"{bn}.gzi" if genome_zip else f"{bn}.fai"
skip_reports = config["skip_reports"]
windowsize  = config["depth_windowsize"]
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
        """
        if (file {input} | grep -q compressed ) ;then
            # is regular gzipped, needs to be BGzipped
            seqtk seq {input} | bgzip -c > {output}
        else
            cp -f {input} {output}
        fi
        """

rule samtools_faidx:
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

rule bwa_index:
    input: 
        f"Genome/{bn}"
    output: 
        multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    log:
        f"Genome/{bn}.bwa.idx.log"
    conda:
        f"{envdir}/align.yaml"
    shell: 
        "bwa index {input} 2> {log}"

rule align:
    input:
        fastq      = get_fq,
        genome     = f"Genome/{bn}",
        genome_idx = multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output:  
        temp(outdir + "/samples/{sample}/{sample}.bwa.sam")
    log:
        outdir + "/logs/bwa/{sample}.bwa.log"
    params: 
        samps = lambda wc: d[wc.get("sample")],
        quality = config["alignment_quality"],
        unmapped = "" if keep_unmapped else "-F 4",
        extra = extra
    benchmark:
        ".Benchmark/Mapping/bwa/align.{sample}.txt"
    threads:
        max(5, workflow.cores - 1)
    conda:
        f"{envdir}/align.yaml"
    shell:
        """
        bwa mem -C -v 2 -t {threads} {params.extra} -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" {input.genome} {input.fastq} 2> {log} |
            samtools view -h {params.unmapped} -q {params.quality} > {output} 
        """

rule mark_duplicates:
    input:
        sam    = outdir + "/samples/{sample}/{sample}.bwa.sam",
        genome = f"Genome/{bn}",
        faidx  = f"Genome/{bn_idx}"
    output:
        temp(outdir + "/samples/{sample}/{sample}.markdup.bam")
    log:
        outdir + "/logs/markdup/{sample}.markdup.log"
    params: 
        tmpdir = lambda wc: outdir + "/." + d[wc.sample]
    resources:
        mem_mb = 2000
    container:
        None
    threads:
        4
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

rule molecule_coverage:
    input:
        stats = outdir + "/reports/data/bxstats/{sample}.bxstats.gz",
        fai = f"Genome/{bn}.fai"
    output: 
        outdir + "/reports/data/coverage/{sample}.molcov.gz"
    params:
        windowsize
    container:
        None
    shell:
        "molecule_coverage.py -f {input.fai} {input.stats} | depth_windows.py {params} | gzip > {output}"

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
        bxstats = outdir + "/reports/data/bxstats/{sample}.bxstats.gz",
        coverage = outdir + "/reports/data/coverage/{sample}.cov.gz",
        molecule_coverage = outdir + "/reports/data/coverage/{sample}.molcov.gz"
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
        outdir + "/reports/bwa.stats.html"
    params:
        outdir = f"{outdir}/reports/data/samtools_stats {outdir}/reports/data/samtools_flagstat",
        options = "--no-version-check --force --quiet --no-data-dir",
        title = "--title \"Basic Alignment Statistics\"",
        comment = "--comment \"This report aggregates samtools stats and samtools flagstats results for all alignments. Samtools stats ignores alignments marked as duplicates.\""
    conda:
        f"{envdir}/qc.yaml"
    shell:
        "multiqc {params} --filename {output} 2> /dev/null"

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
        bams = collect(outdir + "/{sample}.{ext}", sample = samplenames, ext = ["bam", "bam.bai"]),
        reports = collect(outdir + "/reports/{sample}.html", sample = samplenames) if not skip_reports else [],
        agg_report = outdir + "/reports/bwa.stats.html" if not skip_reports else [],
        bx_report = outdir + "/reports/barcodes.summary.html" if (not skip_reports or len(samplenames) == 1) else []
    params:
        quality = config["alignment_quality"],
        unmapped = "" if keep_unmapped else "-F 4",
        extra   = extra
    run:
        summary = ["The harpy align bwa workflow ran using these parameters:"]
        summary.append(f"The provided genome: {genomefile}")
        align = "Sequences were aligned with BWA using:\n"
        align += f"\tbwa mem -C -v 2 {params.extra} -R \"@RG\\tID:SAMPLE\\tSM:SAMPLE\" genome forward_reads reverse_reads |\n"
        align += f"\tsamtools view -h {params.unmapped} -q {params.quality}"
        summary.append(align)
        duplicates = "Duplicates in the alignments were marked following:\n"
        duplicates += "\tsamtools collate |\n"
        duplicates += "\tsamtools fixmate |\n"
        duplicates += f"\tsamtools sort -T SAMPLE --reference {genomefile} -m 2000M |\n"
        duplicates += "\tsamtools markdup -S --barcode-tag BX"
        summary.append(duplicates)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['workflow_call']}"
        summary.append(sm)
        with open(outdir + "/workflow/align.bwa.summary", "w") as f:
            f.write("\n\n".join(summary))