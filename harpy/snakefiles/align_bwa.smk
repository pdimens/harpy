containerized: "docker://pdimens/harpy:latest"

import os
import re
import logging

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake_log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

fqlist       = config["inputs"]["fastq"]
molecule_distance = config["molecule_distance"]
ignore_bx = config["ignore_bx"]
keep_unmapped = config["keep_unmapped"]
extra 		= config.get("extra", "") 
genomefile 	= config["inputs"]["reference"]
bn 			= os.path.basename(genomefile)
workflow_geno = f"workflow/reference/{bn}"
genome_zip  = True if bn.lower().endswith(".gz") else False
workflow_geno_idx = f"{workflow_geno}.gzi" if genome_zip else f"{workflow_geno}.fai"
skip_reports = config["reports"]["skip"]
plot_contigs = config["reports"]["plot_contigs"]    
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
        workflow_geno
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
        workflow_geno
    output: 
        fai = f"{workflow_geno}.fai",
        gzi = f"{workflow_geno}.gzi" if genome_zip else []
    log:
        f"{workflow_geno}.faidx.log"
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

rule make_depth_intervals:
    input:
        fai = f"{workflow_geno}.fai"
    output:
        bed = f"reports/data/coverage/coverage.bed"
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

rule bwa_index:
    input: 
        workflow_geno
    output: 
        multiext(workflow_geno, ".ann", ".bwt", ".pac", ".sa", ".amb")
    log:
        f"{workflow_geno}.bwa.idx.log"
    conda:
        "envs/align.yaml"
    shell: 
        "bwa index {input} 2> {log}"

rule align:
    input:
        fastq      = get_fq,
        genome     = workflow_geno,
        genome_idx = multiext(workflow_geno, ".ann", ".bwt", ".pac", ".sa", ".amb")
    output:
        temp("samples/{sample}/{sample}.bwa.sam")
    log:
        "logs/bwa/{sample}.bwa.log"
    params:
        RG_tag = lambda wc: "\"@RG\\tID:" + wc.get("sample") + "\\tSM:" + wc.get("sample") + "\"",
        samps = lambda wc: d[wc.get("sample")],
        quality = config["alignment_quality"],
        unmapped = "" if keep_unmapped else "-F 4",
        extra = extra
    threads:
        4
    conda:
        "envs/align.yaml"
    shell:
        """
        bwa mem -C -v 2 -t {threads} {params.extra} -R {params.RG_tag} {input.genome} {input.fastq} 2> {log} |
            samtools view -h {params.unmapped} -q {params.quality} > {output} 
        """

rule mark_duplicates:
    input:
        sam    = "samples/{sample}/{sample}.bwa.sam",
        genome = workflow_geno,
        faidx  = workflow_geno_idx
    output:
        temp("samples/{sample}/{sample}.markdup.bam") if not ignore_bx else temp("markdup/{sample}.markdup.bam") 
    log:
        "logs/markdup/{sample}.markdup.log"
    params: 
        tmpdir = lambda wc: "." + d[wc.sample],
        bx_mode = "--barcode-tag BX" if not ignore_bx else ""
    resources:
        mem_mb = 2000
    container:
        None
    threads:
        4
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
        qmd = f"workflow/report/align_stats.qmd"
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
        bam = "{sample}.bam",
        bai = "{sample}.bam.bai"
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
        "reports/bwa.stats.html"
    params:
        outdir = "reports/data/samtools_stats reports/data/samtools_flagstat",
        options = "--no-version-check --force --quiet --no-data-dir",
        title = "--title \"Basic Alignment Statistics\"",
        comment = "--comment \"This report aggregates samtools stats and samtools flagstats results for all alignments. Samtools stats ignores alignments marked as duplicates.\""
    conda:
        "envs/qc.yaml"
    shell:
        "multiqc {params} --filename {output} 2> /dev/null"

rule barcode_report:
    input: 
        f"reports/_quarto.yml",
        f"reports/_harpy.scss",
        collect("reports/data/bxstats/{sample}.bxstats.gz", sample = samplenames),
        qmd = f"workflow/report/align_bxstats.qmd"
    output:
        report = "reports/barcode.summary.html",
        qmd = temp("reports/barcode.summary.qmd")
    params:
        f"reports/data/bxstats/"
    log:
        f"logs/reports/bxstats.report.log"
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
        bams = collect("{sample}.{ext}", sample = samplenames, ext = ["bam", "bam.bai"]),
        samtools = "reports/bwa.stats.html" if not skip_reports else [],
        reports = collect("reports/{sample}.html", sample = samplenames) if not skip_reports and not ignore_bx else [],
        bx_report = "reports/barcode.summary.html" if ((not skip_reports and not ignore_bx) or len(samplenames) == 1) else []
    params:
        quality = config["alignment_quality"],
        unmapped = "" if keep_unmapped else "-F 4",\
        bx_mode = "--barcode-tag BX" if not ignore_bx else "",
        extra   = extra
    run:
        summary = ["The harpy align bwa workflow ran using these parameters:"]
        summary.append(f"The provided genome: {genomefile}")
        align = "Sequences were aligned with BWA using:\n"
        align += f'\tbwa mem -C -v 2 {params.extra} -R "@RG\\tID:SAMPLE\\tSM:SAMPLE" genome forward_reads reverse_reads |\n'
        align += f"\tsamtools view -h {params.unmapped} -q {params.quality}"
        summary.append(align)
        duplicates = "Duplicates in the alignments were marked following:\n"
        duplicates += "\tsamtools collate |\n"
        duplicates += "\tsamtools fixmate |\n"
        duplicates += f"\tsamtools sort -T SAMPLE --reference {genomefile} -m 2000M |\n"
        duplicates += f"\tsamtools markdup -S {params.bx_mode} -d 100 (2500 for novaseq)"
        summary.append(duplicates)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['snakemake_command']}"
        summary.append(sm)
        with open(f"workflow/align.bwa.summary", "w") as f:
            f.write("\n\n".join(summary))