import os

seq_dir = config["seq_directory"]
nbins = config["EMA_bins"]
genomefile = config["genomefile"]
Rsep = config["Rsep"]
fqext = config["fqext"]
samplenames = config["samplenames"]
bed_prox = config["bed_proximity"]
extra = config.get("extra", "") 

bn = os.path.basename(genomefile)
shell("mkdir -p Assembly")
if not os.path.exists(f"Assembly/{bn}"):
    shell(f"ln -sr {genomefile} Assembly/{bn}")

rule create_reports:
    input: 
        expand("Alignments/ema/{sample}.bam", sample = samplenames),
        expand("Alignments/ema/{sample}.bam.bai", sample = samplenames),
        expand("Alignments/ema/stats/moleculesize/{sample}.{ext}", sample = samplenames, ext = ["molsize", "molsize.hist"]),
        expand("Alignments/ema/stats/readsperbx/{sample}.readsperbx", sample = samplenames),
        expand("Alignments/ema/stats/coverage/{sample}.gencov.{aln}.html", sample = samplenames, aln = ["all", "bx"]),
        "Alignments/ema/stats/reads.bxstats.html",
        "Alignments/ema/stats/samtools_stats/alignment.stats.html",
        "Alignments/ema/stats/samtools_flagstat/alignment.flagstat.html"
    message: "Read mapping completed!"
    benchmark: "Benchmark/Mapping/ema/report.txt"
    default_target: True

rule index_genome:
    input: genomefile
    output: 
        asm = f"Assembly/{genomefile}",
        idx = multiext(f"Assembly/{genomefile}", ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
    message: "Indexing {input}"
    log: f"Assembly/{genomefile}.idx.log"
    shell: 
        """
        ln -sr {input} {output.asm}
        bwa index {output.asm} 2> {log}
        samtools faidx --fai-idx {output.asm}.fai {output.asm} 2>> {log}
        """

rule count_beadtags:
    input:
        forward_reads = seq_dir + "/{sample}" + f".{Rsep[0]}.{fqext}",
        reverse_reads = seq_dir + "/{sample}" + f".{Rsep[1]}.{fqext}"
    output: 
        counts = "Alignments/ema/count/{sample}.ema-ncnt",
        logs = temp("Alignments/ema/count/logs/{sample}.count.log")
    wildcard_constraints:
        sample = "[a-zA-Z0-9_-]*"
    message: "Counting barcode frequency: {wildcards.sample}"
    benchmark: "Benchmark/Mapping/ema/Count.{sample}.txt"
    params:
        prefix = lambda wc: "Alignments/ema/count/" + wc.get("sample")
    threads: 1
    shell:
        """
        seqfu interleave -1 {input.forward_reads} -2 {input.reverse_reads} | ema-h count -p -o {params} 2> {output.logs}
        """

rule beadtag_summary:
    input: 
        countlog = expand("Alignments/ema/count/logs/{sample}.count.log", sample = samplenames)
    output: "Alignments/ema/stats/reads.bxstats.html"
    message: "Creating sample barcode validation report"
    benchmark: "Benchmark/Mapping/ema/beadtagsummary.txt"
    script: "../utilities/reportEmaCount.Rmd"

rule preprocess_ema:
    input: 
        forward_reads = seq_dir + "/{sample}" + f".{Rsep[0]}.{fqext}",
        reverse_reads = seq_dir + "/{sample}" + f".{Rsep[1]}.{fqext}",
        emacounts = "Alignments/ema/count/{sample}.ema-ncnt"
    output: 
        bins = temp(expand("Alignments/ema/preproc/{{sample}}/ema-bin-{bin}", bin = ["%03d" % i for i in range(nbins)])),
        unbarcoded = temp("Alignments/ema/preproc/{sample}/ema-nobc")
    wildcard_constraints:
        sample = "[a-zA-Z0-9_-]*"
    log: "Alignments/ema/preproc/logs/{sample}.preproc.log"
    message: "Preprocessing for EMA mapping: {wildcards.sample}"
    benchmark: "Benchmark/Mapping/ema/Preproc.{sample}.txt"
    threads: 2
    params:
        outdir = lambda wc: "Alignments/ema/preproc/" + wc.get("sample"),
        bins = nbins
    shell:
        """
        seqfu interleave -1 {input.forward_reads} -2 {input.reverse_reads} | ema-h preproc -p -n {params.bins} -t {threads} -o {params.outdir} {input.emacounts} 2>&1 | cat - > {log}
        """

rule align_ema:
    input:
        readbin = "Alignments/ema/preproc/{sample}/ema-bin-{bin}",
        genome = f"Assembly/{genomefile}",
        genome_idx = multiext(f"Assembly/{genomefile}", ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
    output: pipe("Alignments/ema/align/{sample}/{sample}.{bin}.sam")
    wildcard_constraints:
        sample = "[a-zA-Z0-9_-]*"
    message: "Mapping onto {input.genome}: {wildcards.sample}-{wildcards.bin}"
    benchmark: "Benchmark/Mapping/ema/Align.{sample}.{bin}.txt"
    params: 
        extra = extra
    threads: 3
    shell:
        """
        ema-h align -t {threads} {params} -d -p haptag -r {input.genome} -o {output} -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" -s {input.readbin} 2> /dev/null
        """

rule align_nobarcode:
    input:
        reads = "Alignments/ema/preproc/{sample}/ema-nobc",
        genome = f"Assembly/{genomefile}",
        genome_idx = multiext(f"Assembly/{genomefile}", ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
    output: 
        samfile = pipe("Alignments/ema/align/{sample}/{sample}.nobarcode.sam")
    benchmark: "Benchmark/Mapping/ema/bwaAlign.{sample}.txt"
    wildcard_constraints:
        sample = "[a-zA-Z0-9_-]*"
    message: "Mapping unbarcoded reads onto {input.genome}: {wildcards.sample}"
    threads: 2
    shell:
        """
        bwa mem -t {threads} -C -M -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" {input.genome} {input.reads} > {output} 2> /dev/null
        """

rule sort_ema:
    input: 
        sam = "Alignments/ema/align/{sample}/{sample}.{emabin}.sam",
        genome = f"Assembly/{genomefile}"
    output: temp("Alignments/ema/align/{sample}/{sample}.{emabin}.bam")
    wildcard_constraints:
        sample = "[a-zA-Z0-9_-]*",
        emabin = "[0-9]*"
    message: "Sorting alignments: {wildcards.sample}-{wildcards.emabin}"
    benchmark: "Benchmark/Mapping/ema/Sort.{sample}.{emabin}.txt"
    threads: 1
    shell: 
        """
        samtools sort -@ {threads} --reference {input.genome} -O bam -l 0 -m 4G -o {output} {input.sam}
        """

rule sort_nobarcode:
    input: 
        sam = "Alignments/ema/align/{sample}/{sample}.nobarcode.sam",
        genome = f"Assembly/{genomefile}"
    output: temp("Alignments/ema/align/{sample}/{sample}.nobarcode.bam.tmp")
    wildcard_constraints:
        sample = "[a-zA-Z0-9_-]*"
    message: "Sorting unbarcoded alignments: {wildcards.sample}"
    benchmark: "Benchmark/Mapping/ema/bwaSort.{sample}.txt"
    threads: 2
    shell:
        """
        samtools sort -@ {threads} -O bam -l 0 -m 4G --reference {input.genome} -o {output} {input.sam}
        """    

rule markduplicates:
    input: "Alignments/ema/align/{sample}/{sample}.nobarcode.bam.tmp"
    output: 
        bam = temp("Alignments/ema/align/{sample}/{sample}.nobarcode.bam"),
        bai = temp("Alignments/ema/align/{sample}/{sample}.nobarcode.bam.bai")
    log: 
        mdlog = "Alignments/ema/stats/markduplicates/{sample}.markdup.nobarcode.log",
        stats = "Alignments/ema/stats/samtools_stats/{sample}.nobarcode.stats",
        flagstat = "Alignments/ema/stats/samtools_flagstat/{sample}.nobarcode.flagstat"
    wildcard_constraints:
        sample = "[a-zA-Z0-9_-]*"
    message: "Marking duplicates in unbarcoded alignments: {wildcards.sample}"
    benchmark: "Benchmark/Mapping/ema/markdup.{sample}.txt"
    threads: 2
    shell:
        """
        sambamba markdup -t {threads} -l 0 {input} {output.bam} 2> {log.mdlog}
        samtools stats {output.bam} > {log.stats}
        samtools flagstat {output.bam} > {log.flagstat}
        """   

rule merge_barcoded:
    input:
        aln_barcoded = expand("Alignments/ema/align/{{sample}}/{{sample}}.{bin}.bam", bin = ["%03d" % i for i in range(nbins)]),
    output: 
        bam = "Alignments/ema/align/barcoded/{sample}.barcoded.bam",
        bai = temp("Alignments/ema/align/barcoded/{sample}.barcoded.bam.bai")
    wildcard_constraints:
        sample = "[a-zA-Z0-9_-]*"
    message: "Merging barcoded alignments: {wildcards.sample}"
    benchmark: "Benchmark/Mapping/ema/merge.{sample}.txt"
    threads: 10
    shell:
        """
        sambamba merge -t {threads} {output.bam} {input}
        """

rule bcstats:
    input: 
        bam = "Alignments/ema/align/barcoded/{sample}.barcoded.bam",
        bai = "Alignments/ema/align/barcoded/{sample}.barcoded.bam.bai"
    log:
        stats = "Alignments/ema/stats/samtools_stats/{sample}.barcoded.stats",
        flagstat = "Alignments/ema/stats/samtools_flagstat/{sample}.barcoded.flagstat"
    wildcard_constraints:
        sample = "[a-zA-Z0-9_-]*"
    message: "Indexing merged barcoded alignemnts: {wildcards.sample}"
    benchmark: "Benchmark/Mapping/ema/indexmerge.{sample}.txt"
    shell:
        """
        samtools stats {input.bam} > {log.stats}
        samtools flagstat {input.bam} > {log.flagstat}
        """

rule BEDconvert:
    input: "Alignments/ema/align/barcoded/{sample}.barcoded.bam"
    output: 
        filt = temp("Alignments/ema/bedfiles/{sample}.bx.bed"),
        unfilt = temp("Alignments/ema/bedfiles/{sample}.all.bed")
    message: "Converting to BED format: {wildcards.sample}"
    wildcard_constraints:
        sample = "[a-zA-Z0-9_-]*"
    params: lambda wc: "Alignments/ema/align/" + wc.get("sample") + "/" + wc.get("sample") + ".barcoded.all.bed"
    threads: 1
    shell:
        """
        utilities/writeBED.pl {input}
        mv {params} {output.unfilt}
        awk '!($4~/A00|B00|C00|D00/)' {output.unfilt} > {output.filt}
        """

rule BX_stats:
    input: "Alignments/ema/bedfiles/{sample}.bx.bed"
    output:	
        molsize = "Alignments/ema/stats/moleculesize/{sample}.molsize",
        molhist = "Alignments/ema/stats/moleculesize/{sample}.molsize.hist",
        readsper = "Alignments/ema/stats/readsperbx/{sample}.readsperbx"
    message: "Calculating molecule size, reads per molecule: {wildcards.sample}"
    wildcard_constraints:
        sample = "[a-zA-Z0-9_-]*"
    threads: 1
    shell:
        """
        cut -f10 {input} | datamash -s groupby 1 count 1 | sort -k 1 -n > {output.readsper}
        awk '{{ print $1"\\t"$2"\\t"$3"\\t"$3-$2"\\t"$4"\\t"$10 }}' {input} | sort -k 4 -n > {output.molsize}
        cut -f4 {output.molsize} | datamash bin:1000 1 | datamash -s groupby 1 count 1 | sort -k 1 -n > {output.molhist}
        # datamash groupby 1,5 min 2 max 3 < {output.molsize} | awk '{{ print $1"\\t"$2"\\t"$4-$3 }}' > 
        """

rule merge_alignments:
    input:
        aln_barcoded = "Alignments/ema/align/barcoded/{sample}.barcoded.bam",
        idx_barcoded = "Alignments/ema/align/barcoded/{sample}.barcoded.bam.bai",
        aln_nobarcode = "Alignments/ema/align/{sample}/{sample}.nobarcode.bam",
        idx_nobarcode = "Alignments/ema/align/{sample}/{sample}.nobarcode.bam.bai"
    output: 
        bam = temp("Alignments/ema/align/{sample}.unsort.bam"),
        bai = temp("Alignments/ema/align/{sample}.unsort.bam.bai")
    wildcard_constraints:
        sample = "[a-zA-Z0-9_-]*"
    message: "Merging all alignments: {wildcards.sample}"
    benchmark: "Benchmark/Mapping/ema/mergebc_nobc.{sample}.txt"
    threads: 10
    shell:
        """
        sambamba merge -t {threads} {output.bam} {input.aln_barcoded} {input.aln_nobarcode}
        """

rule sort_merge:
    input:
        bam = "Alignments/ema/align/{sample}.unsort.bam",
        genome = genomefile
    output: "Alignments/ema/{sample}.bam"
    message: "Sorting merged barcoded alignments: {wildcards.sample}"
    wildcard_constraints:
        sample = "[a-zA-Z0-9_-]*"
    threads: 2
    shell:
        """
        samtools sort -@ {threads} -O bam --reference {input.genome} -l 0 -m 4G -o {output} {input.bam} 2> /dev/null
        """

rule genome_coverage:
    input:
        alntot = "Alignments/ema/{sample}.bam",
        bx = "Alignments/ema/align/barcoded/{sample}.barcoded.bam"
    output: 
        bx = "Alignments/ema/stats/coverage/data/{sample}.bx.gencov",
        alntot = "Alignments/ema/stats/coverage/data/{sample}.all.gencov"
    message: "Calculating genomic coverage: {wildcards.sample}"
    params: bed_prox
    threads: 2
    shell:
        """
        bedtools genomecov -ibam {input.bx}     -bg | bedtools merge -c 4 -o sum -d {params} > {output.bx}
        bedtools genomecov -ibam {input.alntot} -bg | bedtools merge -c 4 -o sum -d {params} > {output.alntot}
        """

rule gencovBX_report:
    input: 
        gencov = "Alignments/ema/stats/coverage/data/{sample}.bx.gencov",
        faidx = f"Assembly/{genomefile}.fai"
    output:
        "Alignments/ema/stats/coverage/{sample}.gencov.bx.html"
    message:
        "Creating report of alignment coverage: {wildcards.sample}"
    script:
        "../utilities/reportGencov.Rmd"

rule gencovAll_report:
    input:
        gencov = "Alignments/ema/stats/coverage/data/{sample}.all.gencov",
        faidx = f"Assembly/{genomefile}.fai"
    output:
        "Alignments/ema/stats/coverage/{sample}.gencov.all.html"
    message:
        "Summarizing alignment coverage: {wildcards.sample}"
    script:
        "../utilities/reportGencov.Rmd"


rule index_alignments:
    input: "Alignments/ema/{sample}.bam"
    output: "Alignments/ema/{sample}.bam.bai"
    message: "Indexing: {input}"
    benchmark: "Benchmark/Mapping/ema/IndexMerged.{sample}.txt"
    wildcard_constraints:
        sample = "[a-zA-Z0-9_-]*"
    shell:
        """
        sambamba index {input} {output}
        """

rule alignment_stats:
    input: 		
        bam = "Alignments/ema/{sample}.bam",
        bai = "Alignments/ema/{sample}.bam.bai"
    output:
        stats = "Alignments/ema/stats/samtools_stats/{sample}.stats",
        flagstat = "Alignments/ema/stats/samtools_flagstat/{sample}.flagstat"
    message: "Calculating alignment stats: {wildcards.sample}"
    benchmark: "Benchmark/Mapping/ema/Mergedstats.{sample}.txt"
    wildcard_constraints:
        sample = "[a-zA-Z0-9_-]*"
    shell:
        """
        samtools stats {input.bam} > {output.stats}
        samtools flagstat {input.bam} > {output.flagstat}
        """

rule samtools_reports:
    input: 
        expand("Alignments/ema/stats/samtools_{ext}/{sample}.{ext}", sample = samplenames, ext = ["stats", "flagstat"]),
    output: 
        stats =    "Alignments/ema/stats/samtools_stats/alignment.stats.html",
        flagstat = "Alignments/ema/stats/samtools_flagstat/alignment.flagstat.html"
    message: "Summarizing samtools stats and flagstats"
    benchmark: "Benchmark/Mapping/ema/report.txt"
    shell:
        """
        multiqc Alignments/ema/stats/samtools_stats    --force --quiet --no-data-dir --filename {output.stats} 2> /dev/null
        multiqc Alignments/ema/stats/samtools_flagstat --force --quiet --no-data-dir --filename {output.flagstat} 2> /dev/null
        """