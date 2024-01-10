import os
import re
import glob
from pathlib import Path

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

conda:
    os.getcwd() + "/harpyenvs/align.yaml"

d = dict(zip(samplenames, samplenames))

def get_fq1(wildcards):
    # code that returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    lst = glob.glob(seq_dir + "/" + wildcards.sample + "*")
    return lst[0]

def get_fq2(wildcards):
    # code that returns a list of fastq files for read 2 based on *wildcards.sample*, e.g.
    lst = sorted(glob.glob(seq_dir + "/" + wildcards.sample + "*"))
    return lst[1]

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
        message:
            "Indexing {input}"
        log:
            f"Genome/{bn}.faidx.gzi.log"
        shell: 
            "samtools faidx --gzi-idx {output.gzi} --fai-idx {output.fai} {input} 2> {log}"
else:
    rule genome_faidx:
        input: 
            f"Genome/{bn}"
        output: 
            f"Genome/{bn}.fai"
        message:
            "Indexing {input}"
        log:
            f"Genome/{bn}.faidx.log"
        shell:
            "samtools faidx --fai-idx {output} {input} 2> {log}"

rule genome_bwa_index:
    input: 
        f"Genome/{bn}"
    output: 
        multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    message:
        "Indexing {input}"
    log:
        f"Genome/{bn}.idx.log"
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

rule beadtag_count:
    input:
        forward_reads = get_fq1,
        reverse_reads = get_fq2
    output: 
        counts = temp(outdir + "/{sample}/{sample}.ema-ncnt"),
        logs   = temp(outdir + "/logs/count/{sample}.count")
    message:
        "Counting barcode frequency: {wildcards.sample}"
    params:
        prefix = lambda wc: outdir + "/" + wc.get("sample") + "/" + wc.get("sample"),
        beadtech = "-p" if platform == "haplotag" else f"-w {whitelist}",
        logdir = f"{outdir}/logs/count/"
    shell:
        """
        mkdir -p {params.prefix} {params.logdir}
        seqtk mergepe {input.forward_reads} {input.reverse_reads} | 
            ema count {params.beadtech} -o {params.prefix} 2> {output.logs}
        """

rule beadtag_summary:
    input: 
        countlog = expand(outdir + "/logs/count/{sample}.count", sample = samplenames)
    output:
        outdir + "/stats/reads.bxcounts.html"
    message:
        "Creating sample barcode validation report"
    script:
        "reportEmaCount.Rmd"

rule preprocess:
    input: 
        forward_reads = get_fq1,
        reverse_reads = get_fq2,
        emacounts     = outdir + "/{sample}/{sample}.ema-ncnt"
    output: 
        bins       	  = temp(expand(outdir + "/{{sample}}/preproc/ema-bin-{bin}", bin = binrange)),
        unbarcoded    = temp(outdir + "/{sample}/preproc/ema-nobc")
    log:
        outdir + "/logs/preproc/{sample}.preproc.log"
    message:
        "Preprocessing for EMA mapping: {wildcards.sample}"
    benchmark:
        ".Benchmark/Mapping/ema/Preproc.{sample}.txt"
    threads:
        2
    params:
        outdir = lambda wc: outdir + "/" + wc.get("sample") + "/preproc",
        beadtech = "-p" if platform == "haplotag" else f"-w {whitelist}",
        bins   = nbins
    shell:
        "seqtk mergepe {input.forward_reads} {input.reverse_reads} | ema preproc {params.beadtech} -n {params.bins} -t {threads} -o {params.outdir} {input.emacounts} 2>&1 | cat - > {log}"

rule align:
    input:
        readbin    = expand(outdir + "/{{sample}}/preproc/ema-bin-{bin}", bin = binrange),
        genome 	   = f"Genome/{bn}",
        geno_faidx = f"Genome/{bn_idx}",
        geno_idx   = multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output:
        aln = temp(outdir + "/{sample}/{sample}.bc.bam"),
        idx = temp(outdir + "/{sample}/{sample}.bc.bam.bai")
    log:
        ema     = outdir + "/logs/{sample}.ema.align.log",
        emasort = outdir + "/logs/{sample}.ema.sort.log"
    message:
        "Aligning barcoded sequences: {wildcards.sample}"
    params: 
        quality = config["quality"],
        tmpdir = lambda wc: outdir + "/." + d[wc.sample],
        beadtech = f"-p {platform}",
        extra = extra
    threads:
        10
    shell:
        """
        EMATHREADS=$(( {threads} - 2 ))
        ema align -t $EMATHREADS {params.extra} -d {params.beadtech} -r {input.genome} -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" -x {input.readbin} 2> {log.ema} |
        samtools view -h -F 4 -q {params.quality} - | 
        samtools sort -T {params.tmpdir} --reference {input.genome} -O bam --write-index -m 4G -o {output.aln}##idx##{output.idx} - 2> {log.emasort}
        rm -rf {params.tmpdir}
        """

rule align_nobarcode:
    input:
        reads      = outdir + "/{sample}/preproc/ema-nobc",
        genome 	   = f"Genome/{bn}",
        geno_faidx = f"Genome/{bn_idx}",
        geno_idx   = multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output: 
        temp(outdir + "/{sample}/{sample}.nobc.bam")
    log:
        bwa     = outdir + "/logs/{sample}.bwa.align.log",
        bwasort = outdir + "/logs/{sample}.bwa.sort.log"
    benchmark:
        ".Benchmark/Mapping/ema/bwaAlign.{sample}.txt"
    params:
        quality = config["quality"]
    message:
        "Aligning unbarcoded sequences: {wildcards.sample}"
    threads:
        8
    shell:
        """
        BWATHREADS=$(( {threads} - 2 ))
        bwa mem -t $BWATHREADS -C -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" {input.genome} {input.reads} 2> {log.bwa} |
        samtools view -h -F 4 -q {params.quality} | 
        samtools sort -O bam -m 4G --reference {input.genome} -o {output} 2> {log.bwasort}
        """

rule mark_duplicates:
    input:
        bam      = outdir + "/{sample}/{sample}.nobc.bam"
    output: 
        bam      = temp(outdir + "/{sample}/{sample}.markdup.nobc.bam"),
        bai      = temp(outdir + "/{sample}/{sample}.markdup.nobc.bam.bai")
    log: 
        mdlog    = outdir + "/logs/markduplicates/{sample}.markdup.nobc.log",
        stats    = outdir + "/stats/samtools_stats/{sample}.nobc.stats",
        flagstat = outdir + "/stats/samtools_flagstat/{sample}.nobc.flagstat"
    message:
        "Marking duplicates in unbarcoded alignments: {wildcards.sample}"
    benchmark:
        ".Benchmark/Mapping/ema/markdup.{sample}.txt"
    threads:
        2
    shell:
        """
        sambamba markdup -t {threads} -l 4 {input} {output.bam} 2> {log.mdlog}
        samtools stats {output.bam} > {log.stats}
        samtools flagstat {output.bam} > {log.flagstat}
        """

rule bx_stats:
    input: 
        bam      = outdir + "/{sample}/{sample}.bc.bam",
        bai      = outdir + "/{sample}/{sample}.bc.bam.bai"
    output:
        stats    = outdir + "/stats/samtools_stats/{sample}.bc.stats",
        flagstat = outdir + "/stats/samtools_flagstat/{sample}.bc.flagstat"
    message:
        "Calculating barcoded alignment stats: {wildcards.sample}"
    shell:
        """
        samtools stats {input.bam} > {output.stats}
        samtools flagstat {input.bam} > {output.flagstat}
        """

rule coverage_stats:
    input: 
        bed     = f"Genome/{bn}.bed",
        nobx    = outdir + "/{sample}/{sample}.markdup.nobc.bam",
        nobxbai = outdir + "/{sample}/{sample}.markdup.nobc.bam.bai",
        bx      = outdir + "/{sample}/{sample}.bc.bam",
        bxbai   = outdir + "/{sample}/{sample}.bc.bam.bai"
    output: 
        outdir + "/stats/coverage/data/{sample}.cov.gz"
    message:
        "Calculating genomic coverage: {wildcards.sample}"
    threads:
        2
    shell:
        "samtools bedcov -c {input.bed} {input.bx} {input.nobx} | gzip > {output}"

rule coverage_report:
    input: 
        outdir + "/stats/coverage/data/{sample}.cov.gz",
    output:
        outdir + "/stats/coverage/{sample}.cov.html"
    message:
        "Creating report of alignment coverage: {wildcards.sample}"
    script:
        "reportEmaGencov.Rmd"

rule merge_alignments:
    input:
        aln_bc   = outdir + "/{sample}/{sample}.bc.bam",
        idx_bc   = outdir + "/{sample}/{sample}.bc.bam.bai",
        aln_nobc = outdir + "/{sample}/{sample}.markdup.nobc.bam",
        idx_nobc = outdir + "/{sample}/{sample}.markdup.nobc.bam.bai"
    output: 
        bam 	 = temp(outdir + "/{sample}/{sample}.unsort.bam"),
        bai 	 = temp(outdir + "/{sample}/{sample}.unsort.bam.bai")
    message:
        "Merging all alignments: {wildcards.sample}"
    benchmark:
        ".Benchmark/Mapping/ema/merge.{sample}.txt"
    threads:
        10
    shell:
        "sambamba merge -t {threads} {output.bam} {input.aln_bc} {input.aln_nobc} 2> /dev/null"

rule sort_merge:
    input:
        bam    = outdir + "/{sample}/{sample}.unsort.bam",
        genome = f"Genome/{bn}"
    output:
        bam = outdir + "/align/{sample}.bam",
        bai = outdir + "/align/{sample}.bam.bai"
    message:
        "Sorting merged barcoded alignments: {wildcards.sample}"
    threads:
        2
    shell:
        "samtools sort -@ {threads} -O bam --reference {input.genome} -m 4G --write-index -o {output.bam}##idx##{output.bai} {input.bam} 2> /dev/null"

rule bx_stats_alignments:
    input:
        bam = outdir + "/align/{sample}.bam",
        bai = outdir + "/align/{sample}.bam.bai"
    output: 
        outdir + "/stats/BXstats/data/{sample}.bxstats.gz"
    message:
        "Calculating barcode alignment statistics: {wildcards.sample}"
    shell:
        "bxStats.py {input.bam} | gzip > {output}"

rule bx_stats_report:
    input:
        outdir + "/stats/BXstats/data/{sample}.bxstats.gz"
    output:	
        outdir + "/stats/BXstats/{sample}.bxstats.html"
    message: 
        "Generating summary of barcode alignment: {wildcards.sample}"
    conda:
        os.getcwd() + "/harpyenvs/r-env.yaml"
    params:
        "none"
    script:
        "reportBxStats.Rmd"

rule general_stats:
    input: 		
        bam      = outdir + "/align/{sample}.bam",
        bai      = outdir + "/align/{sample}.bam.bai"
    output:
        stats    = temp(outdir + "/stats/samtools_stats/{sample}.stats"),
        flagstat = temp(outdir + "/stats/samtools_flagstat/{sample}.flagstat")
    message:
        "Calculating alignment stats: {wildcards.sample}"
    benchmark:
        ".Benchmark/Mapping/ema/Mergedstats.{sample}.txt"
    shell:
        """
        samtools stats {input.bam} > {output.stats}
        samtools flagstat {input.bam} > {output.flagstat}
        """

rule samtools_reports:
    input: 
        expand(outdir + "/stats/samtools_{ext}/{sample}.{ext}", sample = samplenames, ext = ["stats", "flagstat"]),
    output: 
        outdir + "/stats/ema.stats.html"
    message:
        "Summarizing samtools stats and flagstat"
    shell:
        """
        multiqc Align/ema/stats/samtools_stats Align/ema/stats/samtools_flagstat --force --quiet --title "Basic Alignment Statistics" --comment "This report aggregates samtools stats and samtools flagstats results for all alignments." --no-data-dir --filename {output} 2> /dev/null
        """

rule log_runtime:
    output:
        outdir + "/logs/harpy.align.log"
    message:
        "Creating record of relevant runtime parameters: {output}"
    params:
        beadtech = "-p" if platform == "haplotag" else f"-w {whitelist}"
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy align module ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write(f"The directory with sequences: {seq_dir}\n")
            _ = f.write("Barcodes were counted and validated with EMA using:\n")
            _ = f.write("    seqtk mergepe forward.fq.gz reverse.fq.gz | ema count {params.beadtech}\n")
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
            _ = f.write("Overlaps were clipped using:\n")
            _ = f.write("    bam clipOverlap --in file.bam --out outfile.bam --stats --noPhoneHome\n")

rule movelinks:
    default_target: True
    input: 
        bam = expand(outdir + "/align/{sample}.bam", sample = samplenames),
        bai = expand(outdir + "/align/{sample}.bam.bai", sample = samplenames),
        samtools = expand(outdir + "/stats/samtools_{stat}/{sample}.bc.{stat}", stat = ["stats", "flagstat"], sample = samplenames),
        covstats = expand(outdir + "/stats/coverage/{sample}.cov.html", sample = samplenames),
        bxstats = expand(outdir + "/stats/BXstats/{sample}.bxstats.html", sample = samplenames),
        bxcounts = f"{outdir}/stats/reads.bxcounts.html",
        emastats = f"{outdir}/stats/ema.stats.html",
        runlog = f"{outdir}/logs/harpy.align.log"
    message:
        "Finished aligning! Moving alignment files into the base Align/ema directory."
    run:
        for i,j in zip(input.bam, input.bai):
            if not os.path.islink(i):
                # yank out just the filename
                fname = os.path.basename(i)
                # move file into base path
                os.rename(i, f"{outdir}/{fname}")
                # preserve "original" in align folder as symlink
                target = Path(f"{outdir}/{fname}").absolute()
                _ = Path(i).symlink_to(target)
            if not os.path.islink(j):
                # same for .bai file
                fnamebai = os.path.basename(j)
                os.rename(j, f"{outdir}/{fnamebai}")
                targetbai = Path(f"{outdir}/{fnamebai}").absolute()
                _ = Path(j).symlink_to(targetbai)
