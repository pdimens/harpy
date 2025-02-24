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
    sample = r"[a-zA-Z0-9._-]+"

outdir      = config["output_directory"]
nbins 		= config["EMA_bins"]
binrange    = ["%03d" % i for i in range(nbins)]
fqlist       = config["inputs"]["fastq"]
genomefile 	= config["inputs"]["genome"]
platform    = config["platform"]
frag_opt    = config["fragment_density_optimization"]
barcode_list   = config["inputs"].get("barcode_list", "") 
extra 		= config.get("extra", "") 
bn 			= os.path.basename(genomefile)
genome_zip  = True if bn.lower().endswith(".gz") else False
bn_idx      = f"{bn}.gzi" if genome_zip else f"{bn}.fai"
envdir      = os.path.join(os.getcwd(), outdir, "workflow", "envs")
windowsize  = config["depth_windowsize"]
keep_unmapped = config["keep_unmapped"]
skip_reports = config["reports"]["skip"]
plot_contigs = config["reports"]["plot_contigs"]    
bn_r = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
samplenames = {re.sub(bn_r, "", os.path.basename(i), flags = re.IGNORECASE) for i in fqlist}
d = dict(zip(samplenames, samplenames))
os.makedirs(f"{outdir}/logs/ema_count/", exist_ok = True)

def get_fq(wildcards):
    """returns a list of fastq files for read 1 based on *wildcards.sample* e.g."""
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

rule index_genome:
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

rule ema_count:
    input:
        get_fq
    output: 
        counts = temp(outdir + "/ema_count/{sample}.ema-ncnt"),
        logs   = temp(outdir + "/logs/count/{sample}.count")
    params:
        prefix = lambda wc: outdir + "/ema_count/" + wc.get("sample"),
        beadtech = "-p" if platform == "haplotag" else f"-w {barcode_list}"
    conda:
        f"{envdir}/align.yaml"
    shell:
        """
        mkdir -p {params.prefix}
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
        bxtype = "-p" if platform == "haplotag" else f"-w {barcode_list}",
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

rule align_ema:
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
        RG_tag = lambda wc: "\"@RG\\tID:" + wc.get("sample") + "\\tSM:" + wc.get("sample") + "\"",
        bxtype = f"-p {platform}",
        tmpdir = lambda wc: outdir + "/." + d[wc.sample],
        frag_opt = "-d" if frag_opt else "",
        quality = config["alignment_quality"],
        unmapped = "" if keep_unmapped else "-F 4",
        extra = extra
    threads:
        10
    conda:
        f"{envdir}/align.yaml"
    shell:
        """
        ema align -t {threads} {params.extra} {params.frag_opt} {params.bxtype} -r {input.genome} -R {params.RG_tag} -x {input.readbin} 2> {log.ema} |
            samtools view -h {params.unmapped} -q {params.quality} | 
            samtools sort -T {params.tmpdir} --reference {input.genome} -O bam --write-index -m {resources.mem_mb}M -o {output.aln}##idx##{output.idx} - 2> {log.sort}
        rm -rf {params.tmpdir}
        """

rule align_bwa:
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
        quality = config["alignment_quality"],
        unmapped = "" if keep_unmapped else "-F 4",
        RG_tag = lambda wc: "\"@RG\\tID:" + wc.get("sample") + "\\tSM:" + wc.get("sample") + "\""
    threads:
        10
    conda:
        f"{envdir}/align.yaml"
    shell:
        """
        bwa mem -t {threads} -v2 -C -R {params.RG_tag} {input.genome} {input.reads} 2> {log} |
            samtools view -h {params.unmapped} -q {params.quality} > {output}
        """

rule mark_duplicates:
    input:
        sam    = outdir + "/bwa_align/{sample}.bwa.nobc.sam",
        genome = f"Genome/{bn}",
        faidx  = f"Genome/{bn_idx}"
    output:
        temp(outdir + "/bwa_align/{sample}.markdup.nobc.bam")
    log:
        outdir + "/logs/markdup/{sample}.markdup.log"
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
        outdir + "/bwa_align/{sample}.markdup.nobc.bam"
    output:
        temp(outdir + "/bwa_align/{sample}.markdup.nobc.bam.bai")
    container:
        None
    shell:
        "samtools index {input}"

rule concat_alignments:
    priority: 100
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
    shell:
        "samtools depth -a {input.bam} | depth_windows.py {params} | gzip > {output}"

rule barcode_stats:
    input:
        bam = outdir + "/{sample}.bam",
        bai = outdir + "/{sample}.bam.bai"
    output: 
        outdir + "/reports/data/bxstats/{sample}.bxstats.gz"
    container:
        None
    shell:
        "bx_stats.py {input.bam} > {output}"

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
        "molecule_coverage.py -f {input.fai} -w {params} {input.stats} | gzip > {output}"

rule report_config:
    input:
        yaml = f"{outdir}/workflow/report/_quarto.yml",
        scss = f"{outdir}/workflow/report/_harpy.scss"
    output:
        yaml = temp(f"{outdir}/reports/_quarto.yml"),
        scss = temp(f"{outdir}/reports/_harpy.scss")
    run:
        import shutil
        for i,o in zip(input,output):
            shutil.copy(i,o)

rule sample_reports:
    input: 
        f"{outdir}/reports/_quarto.yml",
        f"{outdir}/reports/_harpy.scss",
        bxstats = outdir + "/reports/data/bxstats/{sample}.bxstats.gz",
        coverage = outdir + "/reports/data/coverage/{sample}.cov.gz",
        molecule_coverage = outdir + "/reports/data/coverage/{sample}.molcov.gz",
        qmd = f"{outdir}/workflow/report/align_stats.qmd"
    output:
        report = outdir + "/reports/{sample}.html",
        qmd = temp(outdir + "/reports/{sample}.qmd")
    params:
        mol_dist = f"-P mol_dist:0",
        window_size = f"-P windowsize:{windowsize}",
        contigs = f"-P contigs:{plot_contigs}",
        samplename = lambda wc: "-P sample:" + wc.get("sample")
    log:
        outdir + "/logs/reports/{sample}.alignstats.log"
    conda:
        f"{envdir}/r.yaml"
    shell:
        """
        cp {input.qmd} {output.qmd}
        BXSTATS=$(realpath {input.bxstats})
        COVFILE=$(realpath {input.coverage})
        MOLCOV=$(realpath {input.molecule_coverage})
        quarto render {output.qmd} --log {log} --quiet -P bxstats:$BXSTATS -P coverage:$COVFILE -P molcov:$MOLCOV {params}
        """

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

rule barcode_report:
    input: 
        f"{outdir}/reports/_quarto.yml",
        f"{outdir}/reports/_harpy.scss",
        collect(outdir + "/reports/data/bxstats/{sample}.bxstats.gz", sample = samplenames),
        qmd = f"{outdir}/workflow/report/align_bxstats.qmd"
    output:
        report = f"{outdir}/reports/barcode.summary.html",
        qmd = temp(f"{outdir}/reports/barcode.summary.qmd")
    params:
        f"{outdir}/reports/data/bxstats/"
    log:
        f"{outdir}/logs/reports/bxstats.report.log"
    conda:
        f"{envdir}/r.yaml"
    shell:
        """
        cp {input.qmd} {output.qmd}
        INPATH=$(realpath {params})
        quarto render {output.qmd} --log {log} --quiet -P indir:$INPATH
        """

rule workflow_summary:
    default_target: True
    input:
        bams = collect(outdir + "/{sample}.{ext}", sample = samplenames, ext = [ "bam", "bam.bai"] ),
        cov_report = collect(outdir + "/reports/{sample}.html", sample = samplenames) if not skip_reports else [],
        agg_report = f"{outdir}/reports/ema.stats.html" if not skip_reports else [],
        bx_report = outdir + "/reports/barcode.summary.html" if (not skip_reports or len(samplenames) == 1) else []
    params:
        beadtech = "-p" if platform == "haplotag" else f"-w {barcode_list}",
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
        ema_align += f'\tema align {extra} {params.frag_opt} -p {platform} -R "@RG\\tID:SAMPLE\\tSM:SAMPLE" |\n'
        ema_align += f"\tsamtools view -h {params.unmapped} -q {config["alignment_quality"]} - |\n"
        ema_align += "\tsamtools sort --reference genome"
        summary.append(ema_align)
        bwa_align = "Non-barcoded and invalid-barcoded sequences were aligned with BWA using:\n"
        bwa_align += '\tbwa mem -C -v 2 -R "@RG\\tID:SAMPLE\\tSM:SAMPLE" genome forward_reads reverse_reads |\n'
        bwa_align += f"\tsamtools view -h {params.unmapped} -q {config["alignment_quality"]}"
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
        sm += f"\t{config['workflow_call']}"
        summary.append(sm)
        with open(outdir + "/workflow/align.ema.summary", "w") as f:
            f.write("\n\n".join(summary))
