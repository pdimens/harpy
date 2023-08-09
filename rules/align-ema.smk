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
extra 		= config.get("extra", "") 
bn 			= os.path.basename(genomefile)

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
        "Symlinking {input} to Genome/"
    shell: 
        "ln -sr {input} {output}"

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
        counts = outdir + "/{sample}/{sample}.ema-ncnt",
        logs   = temp(outdir + "/logs/count/{sample}.count")
    message:
        "Counting barcode frequency: {wildcards.sample}"
    params:
        prefix = lambda wc: outdir + "/" + wc.get("sample") + "/" + wc.get("sample")
    shell:
        "seqfu interleave -1 {input.forward_reads} -2 {input.reverse_reads} | ema count -p -o {params} 2> {output.logs}"

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
        "Benchmark/Mapping/ema/Preproc.{sample}.txt"
    threads:
        2
    params:
        outdir = lambda wc: outdir + "/" + wc.get("sample") + "/preproc",
        bins   = nbins
    shell:
        "seqfu interleave -1 {input.forward_reads} -2 {input.reverse_reads} | ema preproc -p -n {params.bins} -t {threads} -o {params.outdir} {input.emacounts} 2>&1 | cat - > {log}"


rule align:
    input:
        readbin  = expand(outdir + "/{{sample}}/preproc/ema-bin-{bin}", bin = binrange),
        genome 	 = f"Genome/{bn}",
        geno_idx = multiext(f"Genome/{bn}", ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
    output:
        aln = temp(outdir + "/{sample}/{sample}.bc.bam"),
        idx = temp(outdir + "/{sample}/{sample}.bc.bam.bai")
    message:
        "Aligning barcoded sequences: {wildcards.sample}"
    params: 
        quality = config["quality"],
        extra = extra
    threads:
        8
    shell:
        """
        EMATHREADS=$(( {threads} - 2 ))
        ema align -t $EMATHREADS {params.extra} -d -p haplotag -r {input.genome} -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" -x {input.readbin} 2> /dev/null |
        samtools view -h -F 4 -q {params.quality} - | 
        samtools sort --reference {input.genome} -O bam --write-index -m 4G -o {output.aln}##idx##{output.idx} - 2> /dev/null
        """

rule align_nobarcode:
    input:
        reads      = outdir + "/{sample}/preproc/ema-nobc",
        genome 	   = f"Genome/{bn}",
        genome_idx = multiext(f"Genome/{bn}", ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
    output: 
        temp(outdir + "/{sample}/{sample}.nobc.bam")
    benchmark:
        "Benchmark/Mapping/ema/bwaAlign.{sample}.txt"
    params:
        quality = config["quality"]
    message:
        "Aligning unbarcoded sequences: {wildcards.sample}"
    threads:
        8
    shell:
        """
        BWATHREADS=$(( {threads} - 2 ))
        bwa mem -t $BWATHREADS -C -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" {input.genome} {input.reads} 2> /dev/null |
        samtools view -h -F 4 -q {params.quality} | 
        samtools sort -O bam -m 4G --reference {input.genome} -o {output} 2> /dev/null
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
        "Benchmark/Mapping/ema/markdup.{sample}.txt"
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
        "Benchmark/Mapping/ema/merge.{sample}.txt"
    threads:
        10
    shell:
        "sambamba merge -t {threads} {output.bam} {input.aln_bc} {input.aln_nobc} 2> /dev/null"

rule sort_merge:
    input:
        bam    = outdir + "/{sample}/{sample}.unsort.bam",
        genome = f"Genome/{bn}"
    output:
        bam = temp(outdir + "/{sample}/{sample}.sorted.bam"),
        bai = temp(outdir + "/{sample}/{sample}.sorted.bam.bai")
    message:
        "Sorting merged barcoded alignments: {wildcards.sample}"
    threads:
        2
    shell:
        "samtools sort -@ {threads} -O bam --reference {input.genome} -m 4G --write-index -o {output.bam}##idx##{output.bai} {input.bam} 2> /dev/null"

rule bx_stats_alignments:
    input:
        bam = outdir + "/{sample}/{sample}.sorted.bam",
        bai = outdir + "/{sample}/{sample}.sorted.bam.bai"
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
    script:
        "reportBxStats.Rmd"

rule clip_overlap:
    input:
        bam = outdir + "/{sample}/{sample}.sorted.bam",
        bai = outdir + "/{sample}/{sample}.sorted.bam.bai"
    output:
        bam = outdir + "/align/{sample}.bam",
        bai = outdir + "/align/{sample}.bam.bai"
    log:
        outdir + "/logs/clipOverlap/{sample}.clipOverlap.log"
    message:
        "Clipping alignment overlaps: {wildcards.sample}"
    shell:
        """
        bam clipOverlap --in {input.bam} --out {output.bam} --stats --noPhoneHome > {log} 2>&1
        samtools index {output.bam}
        """

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
        "Benchmark/Mapping/ema/Mergedstats.{sample}.txt"
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
        "multiqc Align/ema/stats/samtools_stats Align/ema/stats/samtools_flagstat --force --quiet --no-data-dir --filename {output} 2> /dev/null"

rule log_runtime:
    output:
        outdir + "/logs/harpy.align.log"
    message:
        "Creating record of relevant runtime parameters: {output}"
    params:
        extra = extra
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy align module ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write(f"The directory with sequences: {seq_dir}\n")
            _ = f.write("Barcodes were counted and validated with EMA using:\n")
            _ = f.write("\tseqfu interleave forward.fq.gz reverse.fq.gz | ema count -p\n")
            _ = f.write("Barcoded sequences were binned with EMA using:\n")
            _ = f.write(f"\tseqfu interleave forward.fq.gz reverse.fq.gz | ema preproc -p -n {nbins}\n")
            _ = f.write("Barcoded bins were aligned with ema align using:\n")
            _ = f.write("\tema align " + extra + " -d -p haptag -R \"@RG\\tID:SAMPLE\\tSM:SAMPLE\" |\n")
            _ = f.write("\tsamtools view -h -F 4 -q " + str(config["quality"]) + " - |\n") 
            _ = f.write("\tsamtools sort --reference genome -m 4G\n\n")
            _ = f.write("Invalid/non barcoded sequences were aligned with BWA using:\n")
            _ = f.write("\tbwa mem -C -R \"@RG\\tID:SAMPLE\\tSM:SAMPLE\" genome forward_reads reverse_reads |\n")
            _ = f.write("\tsambamba markdup -l 0\n")

rule all:
    input: 
        bam = expand(outdir + "/align/{sample}.bam", sample = samplenames),
        samtools = expand(outdir + "/stats/samtools_{stat}/{sample}.bc.{stat}", stat = ["stats", "flagstat"], sample = samplenames),
        covstats = expand(outdir + "/stats/coverage/{sample}.cov.html", sample = samplenames),
        bxstats = expand(outdir + "/stats/BXstats/{sample}.bxstats.html", sample = samplenames),
        bxcounts = f"{outdir}/stats/reads.bxcounts.html",
        emastats = f"{outdir}/stats/ema.stats.html",
        runlog = f"{outdir}/logs/harpy.align.log"
        #expand(outdir + "/{sample}/{sample}.markdup.nobc.bam", sample = samplenames),
        #expand(outdir + "/{sample}/{sample}.bc.bam", sample = samplenames)
    message:
        "Finished aligning! Moving alignment files into the base Align/ema directory."
    default_target: True
    run:
        for i in input[0]:
            fname = os.path.basename(i)
            try:
                # move file into base path
                os.rename(i, f"{outdir}/{fname}")
                # preserve "original" in align folder as symlink
                target = Path(f"{outdir}/{fname}").absolute()
                _ = Path(i).symlink_to(target)
            except:
                pass
