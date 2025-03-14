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
workflowdir = f"{workflowdir}"
envdir      = os.path.join(os.getcwd(), outdir, "workflow", "envs")
fqlist      = config["inputs"]["fastq"]
extra 		= config.get("extra", "") 
genomefile 	= config["inputs"]["genome"]
bn 			= os.path.basename(genomefile)
if bn.lower().endswith(".gz"):
    bn = bn[:-3]
workflow_geno = f"{workflowdir}/genome/{bn}"
windowsize  = config["depth_windowsize"]
molecule_distance = config["molecule_distance"]
ignore_bx = config["ignore_bx"]
keep_unmapped = config["keep_unmapped"]
readlen = config["average_read_length"]
autolen = isinstance(readlen, str)
skip_reports = config["reports"]["skip"]
plot_contigs = config["reports"]["plot_contigs"]    
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
        "seqtk seq {input} > {output}"

rule index_genome:
    input: 
        workflow_geno
    output: 
        f"{workflow_geno}.fai",
    log:
        f"{workflow_geno}.faidx.log"
    container:
        None
    shell: 
        "samtools faidx --fai-idx {output} {input} 2> {log}"

rule make_depth_intervals:
    input:
        f"{workflow_geno}.fai"
    output:
        outdir + "/reports/data/coverage/coverage.bed"
    run:
        with open(input[0], "r") as fai, open(output[0], "w") as bed:
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

rule strobe_index:
    input: 
        workflow_geno
    output:
        f"{workflow_geno}.r{readlen}.sti"
    log:
        f"{workflow_geno}.r{readlen}.sti.log"
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
        genome   = workflow_geno,
        genome_index   = f"{workflow_geno}.r{readlen}.sti" if not autolen else []
    output:  
        temp(outdir + "/samples/{sample}/{sample}.strobe.sam")
    log:
        outdir + "/logs/strobealign/{sample}.strobealign.log"
    params: 
        samps = lambda wc: d[wc.get("sample")],
        readlen = "" if autolen else f"--use-index -r {readlen}",
        quality = config["alignment_quality"],
        unmapped_strobe = "" if keep_unmapped else "-U",
        unmapped = "" if keep_unmapped else "-F 4",
        extra = extra
    threads:
        4
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
        genome = workflow_geno,
        faidx  = f"{workflow_geno}.fai"
    output:
        temp(outdir + "/samples/{sample}/{sample}.markdup.bam") if not ignore_bx else temp(outdir + "/markdup/{sample}.markdup.bam")
    log:
        outdir + "/logs/markdup/{sample}.markdup.log"
    params: 
        tmpdir = lambda wc: outdir + "/." + d[wc.sample],
        bx_mode = "--barcode-tag BX" if not ignore_bx else ""
    resources:
        mem_mb = 2000
    threads:
        2
    container:
        None
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
        bam = outdir + "/samples/{sample}/{sample}.markdup.bam",
    output:
        bam = outdir + "/{sample}.bam",
        bai = outdir + "/{sample}.bam.bai"
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
        bam = outdir + "/{sample}.bam",
        bai = outdir + "/{sample}.bam.bai"
    output: 
        outdir + "/reports/data/bxstats/{sample}.bxstats.gz"
    params:
        sample = lambda wc: d[wc.sample]
    container:
        None
    shell:
        "bx_stats.py {input.bam} > {output}"

rule molecule_coverage:
    input:
        stats = outdir + "/reports/data/bxstats/{sample}.bxstats.gz",
        fai = f"{workflow_geno}.fai"
    output: 
        outdir + "/reports/data/coverage/{sample}.molcov.gz"
    params:
        windowsize
    container:
        None
    shell:
        "molecule_coverage.py -f {input.fai} -w {params} {input.stats} | gzip > {output}"

rule alignment_coverage:
    input: 
        bam = outdir + "/{sample}.bam",
        bai = outdir + "/{sample}.bam.bai",
        bed = outdir + "/reports/data/coverage/coverage.bed"
    output: 
        outdir + "/reports/data/coverage/{sample}.cov.gz"
    container:
        None
    shell:
        "samtools bedcov -c {input.bed} {input.bam} | awk '{{ $5 = $5 / ($3 + 1 - $2); print }}' | gzip > {output}"

rule report_config:
    input:
        yaml = f"{workflowdir}/report/_quarto.yml",
        scss = f"{workflowdir}/report/_harpy.scss"
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
        qmd = f"{workflowdir}/report/align_stats.qmd"
    output:
        report = outdir + "/reports/{sample}.html",
        qmd = temp(outdir + "/reports/{sample}.qmd")
    params:
        mol_dist = f"-P mol_dist:{molecule_distance}",
        window_size = f"-P windowsize:{windowsize}",
        contigs = f"-P contigs:{plot_contigs}",
        samplename = lambda wc: "-P sample:" + wc.get("sample")
    log:
        outdir + "/logs/reports/{sample}.alignstats.log"
    conda:
        f"{envdir}/r.yaml"
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
            bam = outdir + "/markdup/{sample}.markdup.bam"
        output:
            bam = outdir + "/{sample}.bam",
            bai = outdir + "/{sample}.bam.bai"
        container:
            None
        shell:
            """
            mv {input.bam} {output.bam}
            samtools index {output.bam}
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
        f"{outdir}/reports/_quarto.yml",
        f"{outdir}/reports/_harpy.scss",
        collect(outdir + "/reports/data/bxstats/{sample}.bxstats.gz", sample = samplenames),
        qmd = f"{workflowdir}/report/align_bxstats.qmd"
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
        cp -f {input.qmd} {output.qmd}
        INPATH=$(realpath {params})
        quarto render {output.qmd} --log {log} --quiet -P indir:$INPATH
        """

rule workflow_summary:
    default_target: True
    input: 
        bams = collect(outdir + "/{sample}.{ext}", sample = samplenames, ext = ["bam","bam.bai"]),
        samtools =  outdir + "/reports/strobealign.stats.html" if not skip_reports else [] ,
        reports = collect(outdir + "/reports/{sample}.html", sample = samplenames) if not skip_reports and not ignore_bx else [],
        bx_report = outdir + "/reports/barcode.summary.html" if ((not skip_reports and not ignore_bx) or len(samplenames) == 1) else []
    params:
        readlen = readlen,
        quality = config["alignment_quality"],
        unmapped_strobe = "" if keep_unmapped else "-U",
        unmapped = "" if keep_unmapped else "-F 4",
        bx_mode = "--barcode-tag BX" if not ignore_bx else "",
        extra   = extra
    run:
        summary = ["The harpy align strobe workflow ran using these parameters:"]
        summary.append(f"The provided genome: {genomefile}")
        align = "Sequences were aligned with strobealign using:\n"
        if autolen:
            align += f"\tstrobealign -U -C --rg-id=SAMPLE --rg=SM:SAMPLE {params.extra} genome reads.F.fq reads.R.fq |\n"
        else:
            align = "The genome index was created using:\n"
            align += f"\tstrobealign --create-index -r {params.readlen} genome\n\n"
            align += "Sequences were aligned with strobealign using:\n"
            align += f"\tstrobealign --use-index {params.unmapped_strobe} -N 2 -C --rg=SM:SAMPLE {params.extra} genome reads.F.fq reads.R.fq |\n"
        align += f"\t\tsamtools view -h {params.unmapped} -q {params.quality}"
        summary.append(align)
        duplicates = "Duplicates in the alignments were marked following:\n"
        duplicates += "\tsamtools collate |\n"
        duplicates += "\tsamtools fixmate |\n"
        duplicates += f"\tsamtools sort -T SAMPLE --reference {genomefile} -m 2000M |\n"
        duplicates += f"\tsamtools markdup -S {params.bx_mode} -d 100 (2500 for novaseq)"
        summary.append(duplicates)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['workflow_call']}"
        summary.append(sm)
        with open(f"{workflowdir}/align.strobealign.summary", "w") as f:
            f.write("\n\n".join(summary))