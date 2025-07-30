containerized: "docker://pdimens/harpy:latest"

import os
import re
import logging

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake"]["log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

nbins 		= config["EMA_bins"]
binrange    = ["%03d" % i for i in range(nbins)]
fqlist       = config["inputs"]["fastq"]
lr_platform    = config["platform"]
lr_platform = "haplotag" if lr_platform == "haplotagging" else lr_platform
frag_opt    = config["fragment_density_optimization"]
barcode_list   = config["inputs"].get("barcode_list", "") 
extra 		= config.get("extra", "") 
genomefile 	= config["inputs"]["reference"]
bn 			= os.path.basename(genomefile)
workflow_geno = f"workflow/reference/{bn}"
genome_zip  = True if bn.lower().endswith(".gz") else False
workflow_geno_idx = f"{workflow_geno}.gzi" if genome_zip else f"{workflow_geno}.fai"
windowsize  = config["depth_windowsize"]
keep_unmapped = config["keep_unmapped"]
skip_reports = config["reports"]["skip"]
plot_contigs = config["reports"]["plot_contigs"]    
bn_r = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
samplenames = {re.sub(bn_r, "", os.path.basename(i), flags = re.IGNORECASE) for i in fqlist}
d = dict(zip(samplenames, samplenames))
os.makedirs("logs/ema_count/", exist_ok = True)

def get_fq(wildcards):
    """returns a list of fastq files for read 1 based on *wildcards.sample* e.g."""
    r = re.compile(fr".*/({re.escape(wildcards.sample)}){bn_r}", flags = re.IGNORECASE)
    return sorted(list(filter(r.match, fqlist))[:2])

rule preprocess_reference:
    input:
        genomefile
    output: 
        geno = workflow_geno,
        bwa_idx = multiext(workflow_geno, ".ann", ".bwt", ".pac", ".sa", ".amb"),
        fai = f"{workflow_geno}.fai",
        gzi = f"{workflow_geno}.gzi" if genome_zip else []
    log:
        f"{workflow_geno}.preprocess.log"
    params:
        genome_zip
    conda:
        "envs/align.yaml"
    shell: 
        """
        if (file {input} | grep -q compressed ) ;then
            # is regular gzipped, needs to be BGzipped
            seqtk seq {input} | bgzip -c > {output.geno}
        else
            cp -f {input} {output.geno}
        fi

        if [ "{params}" = "True" ]; then
            samtools faidx --gzi-idx {output.gzi} --fai-idx {output.fai} {output.geno} 2>> {log}
        else
            samtools faidx --fai-idx {output.fai} {output.geno} 2>> {log}
        fi

        bwa index {output.geno} 2> {log}
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

rule ema_count:
    input:
        get_fq
    output: 
        temp("ema_count/{sample}.ema-ncnt"),
        logs   = temp("logs/count/{sample}.count")
    params:
        prefix = lambda wc: "ema_count/" + wc.get("sample"),
        beadtech = "-p" if lr_platform == "haplotag" else f"-w {barcode_list}"
    conda:
        "envs/align.yaml"
    shell:
        """
        mkdir -p {params.prefix}
        seqtk mergepe {input} |
        ema count {params.beadtech} -o {params.prefix} 2> {output.logs}
        """

rule ema_preprocess:
    input: 
        reads = get_fq,
        emacounts  = "ema_count/{sample}.ema-ncnt"
    output: 
        temp("ema_preproc/{sample}/ema-nobc"),
        bins = temp(collect("ema_preproc/{{sample}}/ema-bin-{bin}", bin = binrange))
    log:
        "logs/ema_preproc/{sample}.preproc.log"
    params:
        outdir = lambda wc: "ema_preproc/" + wc.get("sample"),
        bxtype = "-p" if lr_platform == "haplotag" else f"-w {barcode_list}",
        bins   = nbins
    threads:
        2
    conda:
        "envs/align.yaml"
    shell:
        """
        seqtk mergepe {input.reads} |
        ema preproc {params.bxtype} -n {params.bins} -t {threads} -o {params.outdir} {input.emacounts} 2>&1 |
        cat - > {log}
        """

rule align_ema:
    input:
        readbin    = collect("ema_preproc/{{sample}}/ema-bin-{bin}", bin = binrange),
        genome 	   = workflow_geno,
        geno_faidx = workflow_geno_idx,
        geno_idx   = multiext(workflow_geno, ".ann", ".bwt", ".pac", ".sa", ".amb")
    output:
        aln = temp("ema_align/{sample}.ema.bam")
    log:
        ema  = "logs/align/{sample}.ema.align.log",
        sort = "logs/align/{sample}.ema.sort.log"
    resources:
        mem_mb = 500
    params:
        RG_tag = lambda wc: "\"@RG\\tID:" + wc.get("sample") + "\\tSM:" + wc.get("sample") + "\"",
        bxtype = f"-p {lr_platform}",
        tmpdir = lambda wc: "." + d[wc.sample],
        frag_opt = "-d" if frag_opt else "",
        quality = config["alignment_quality"],
        unmapped = "" if keep_unmapped else "-F 4",
        extra = extra
    threads:
        10
    conda:
        "envs/align.yaml"
    shell:
        """
        ema align -t {threads} {params.extra} {params.frag_opt} {params.bxtype} -r {input.genome} -R {params.RG_tag} -x {input.readbin} 2> {log.ema} |
            samtools view -h {params.unmapped} -q {params.quality} | 
            samtools sort -T {params.tmpdir} --reference {input.genome} -O bam -m {resources.mem_mb}M -o {output.aln} - 2> {log.sort}
        rm -rf {params.tmpdir}
        """

rule standardize_barcodes:
    input:
        "ema_align/{sample}.ema.bam"
    output:
        temp("ema_align/{sample}.bc.bam.bai"),
        bam = temp("ema_align/{sample}.bc.bam")
    log:
        "logs/{sample}.standardize.log"
    container:
        None
    shell:
        """
        standardize_barcodes_sam.py < {input} 2> {log} | samtools view -h -b > {output.bam}
        samtools index {output.bam}
        """

rule align_bwa:
    input:
        reads      = "ema_preproc/{sample}/ema-nobc",
        genome 	   = workflow_geno,
        geno_faidx = workflow_geno_idx,
        geno_idx   = multiext(workflow_geno, ".ann", ".bwt", ".pac", ".sa", ".amb")
    output: 
        temp("bwa_align/{sample}.bwa.nobc.sam")
    log:
        "logs/align/{sample}.bwa.align.log"
    params:
        quality = f"-T {config['alignment_quality']}",
        unmapped = "" if keep_unmapped else "| samtools view -h -F 4",
        RG_tag = lambda wc: "-R \"@RG\\tID:" + wc.get("sample") + "\\tSM:" + wc.get("sample") + "\""
    threads:
        10
    conda:
        "envs/align.yaml"
    shell:
        """
        bwa mem -t {threads} -v2 -C {params.quality} {params.RG_tag} {input.genome} {input.reads} 2> {log} {params.unmapped} > {output}
        """

rule mark_duplicates:
    input:
        sam    = "bwa_align/{sample}.bwa.nobc.sam",
        genome = workflow_geno,
        faidx  = workflow_geno_idx
    output:
        temp("bwa_align/{sample}.markdup.nobc.bam")
    log:
        "logs/markdup/{sample}.markdup.log"
    params: 
        tmpdir = lambda wc: "." + d[wc.sample]
    resources:
        mem_mb = 500
    container:
        None
    threads:
        2
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
            samtools markdup -@ {threads} -S -d $OPTICAL_BUFFER -f {log} - {output}
        rm -rf {params.tmpdir}
        """

rule index_duplicates:
    input:
        "bwa_align/{sample}.markdup.nobc.bam"
    output:
        temp("bwa_align/{sample}.markdup.nobc.bam.bai")
    container:
        None
    shell:
        "samtools index {input}"

rule concat_alignments:
    priority: 100
    input:
        aln_bc   = "ema_align/{sample}.bc.bam",
        idx_bc   = "ema_align/{sample}.bc.bam.bai",
        aln_nobc = "bwa_align/{sample}.markdup.nobc.bam",
        idx_nobc = "bwa_align/{sample}.markdup.nobc.bam.bai",
        genome   = workflow_geno
    output: 
        bam = "{sample}.bam",
        bai = "{sample}.bam.bai"
    log:
        "logs/concat/{sample}.concat.log"
    resources:
        mem_mb = 500
    container:
        None
    shell:
        """
        samtools cat {input.aln_bc} {input.aln_nobc} |
            samtools sort -@ 1 -O bam --reference {input.genome} -m {resources.mem_mb}M --write-index -o {output.bam}##idx##{output.bai} - 2> {log}
        """

rule alignment_coverage:
    input: 
        "{sample}.bam.bai",
        bam = "{sample}.bam",
        bed = "reports/data/coverage/coverage.bed"
    output: 
        "reports/data/coverage/{sample}.cov.gz"
    container:
        None
    shell:
        "samtools bedcov -c {input.bed} {input.bam} | awk '{{ $6 = ($4 / ($3 + 1 - $2)); print }}' | gzip > {output}"

rule barcode_stats:
    input:
        "{sample}.bam.bai",
        bam = "{sample}.bam"
    output: 
        "reports/data/bxstats/{sample}.bxstats.gz"
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
    log:
        "logs/molcov/{sample}.molcov.log"
    params:
        windowsize
    container:
        None
    shell:
        "molecule_coverage.py -f {input.fai} -w {params} {input.stats} 2> {log} | gzip > {output}"

rule configure_report:
    input:
        yaml = "workflow/report/_quarto.yml",
        scss = "workflow/report/_harpy.scss"
    output:
        yaml = temp(f"reports/_quarto.yml"),
        scss = temp(f"reports/_harpy.scss")
    run:
        import shutil
        for i,o in zip(input,output):
            shutil.copy(i,o)

rule sample_reports:
    input: 
        f"reports/_quarto.yml",
        f"reports/_harpy.scss",
        bxstats = "reports/data/bxstats/{sample}.bxstats.gz",
        coverage = "reports/data/coverage/{sample}.cov.gz",
        molecule_coverage = "reports/data/coverage/{sample}.molcov.gz",
        qmd = "workflow/report/align_stats.qmd"
    output:
        report = "reports/{sample}.html",
        qmd = temp("reports/{sample}.qmd")
    params:
        mol_dist = f"-P mol_dist:0",
        window_size = f"-P windowsize:{windowsize}",
        contigs = f"-P contigs:{plot_contigs}",
        samplename = lambda wc: "-P sample:" + wc.get("sample")
    log:
        "logs/reports/{sample}.alignstats.log"
    retries:
        3
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

rule general_stats:
    input: 		
        "{sample}.bam.bai",
        bam = "{sample}.bam"
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
        collect("reports/data/samtools_{ext}/{sample}.{ext}", sample = samplenames, ext = ["stats", "flagstat"]),
    output: 
        "reports/ema.stats.html"
    params:
        outdir = f"reports/data/samtools_stats reports/data/samtools_flagstat",
        options = "--no-ai --no-version-check --force --quiet --no-data-dir",
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
        qmd = "workflow/report/align_bxstats.qmd"
    output:
        report = f"reports/barcode.summary.html",
        qmd = temp(f"reports/barcode.summary.qmd")
    params:
        f"reports/data/bxstats/"
    log:
        f"logs/reports/bxstats.report.log"
    conda:
        "envs/r.yaml"
    retries:
        3
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        INPATH=$(realpath {params})
        quarto render {output.qmd} --log {log} --quiet -P indir:$INPATH
        """

rule workflow_summary:
    default_target: True
    input:
        bams = collect("{sample}.{ext}", sample = samplenames, ext = [ "bam", "bam.bai"] ),
        cov_report = collect("reports/{sample}.html", sample = samplenames) if not skip_reports else [],
        agg_report = f"reports/ema.stats.html" if not skip_reports else [],
        bx_report = "reports/barcode.summary.html" if (not skip_reports and len(samplenames) > 1) else []
    params:
        beadtech = "-p" if lr_platform == "haplotag" else f"-w {barcode_list}",
        unmapped = "" if keep_unmapped else "-F 4",
        frag_opt = "-d" if frag_opt else ""
    run:
        summary = ["The harpy align ema workflow ran using these parameters:"]
        summary.append(f"The provided genome: {genomefile}")
        counts = "Barcodes were counted and validated with EMA using:\n"
        counts += f"\tseqtk mergepe forward.fq.gz reverse.fq.gz | ema count {params.beadtech}"
        summary.append(counts)
        bins = "Barcoded sequences were binned with EMA using:\n"
        bins += f"\tseqtk mergepe forward.fq.gz reverse.fq.gz | ema preproc {params.beadtech} -n {nbins}"
        summary.append(bins)
        ema_align = "Barcoded bins were aligned with ema align using:\n"
        ema_align += f'\tema align {extra} {params.frag_opt} -p {lr_platform} -R "@RG\\tID:SAMPLE\\tSM:SAMPLE" |\n'
        ema_align += f"\tsamtools view -h {params.unmapped} -q {config["alignment_quality"]} - |\n"
        ema_align += "\tsamtools sort --reference genome"
        summary.append(ema_align)
        bwa_align = "Non-barcoded and invalid-barcoded sequences were aligned with BWA using:\n"
        bwa_align += '\tbwa mem -C -v 2 -R "@RG\\tID:SAMPLE\\tSM:SAMPLE" genome forward_reads reverse_reads |\n'
        bwa_align += f"\tsamtools view -h {params.unmapped} -q {config["alignment_quality"]}"
        standardization = "Barcodes were standardized in the EMA aligments using:\n"
        standardization += "\tstandardize_barcodes_sam.py < {input} | samtools view -h -b > {output}"
        summary.append(standardization)
        summary.append(bwa_align)
        duplicates = "Duplicates in non-barcoded alignments were marked following:\n"
        duplicates += "\tsamtools collate |\n"
        duplicates += "\tsamtools fixmate |\n"
        duplicates += f"\tsamtools sort -T SAMPLE --reference {genomefile} |\n"
        duplicates += "\tsamtools markdup -S -d 100 (2500 for novaseq)"
        summary.append(duplicates)
        merged = "Alignments were merged using:\n"
        merged += "\tsamtools cat barcode.bam nobarcode.bam > concat.bam"
        summary.append(merged)
        sorting = "Merged alignments were sorted using:\n"
        sorting += "\tsamtools sort -m 2000M concat.bam"
        summary.append(sorting)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['snakemake']['relative']}"
        summary.append(sm)
        with open("workflow/align.ema.summary", "w") as f:
            f.write("\n\n".join(summary))
