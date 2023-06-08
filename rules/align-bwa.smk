import os

seq_dir		= config["seq_directory"]
genomefile 	= config["genomefile"]
Rsep 		= config["Rsep"]
fqext 		= config["fqext"]
samplenames = config["samplenames"]
extra 		= config.get("extra", "") 
bn 			= os.path.basename(genomefile)
outdir      = "Align/bwa"

rule create_reports:
	input: 
		expand(outdir + "/{sample}.bam", sample = samplenames),
		expand(outdir + "/stats/coverage/{sample}.cov.html", sample = samplenames),
		expand(outdir + "/stats/moleculesize/{sample}.{ext}", sample = samplenames, ext = ["molsize.gz", "molsize.hist"]),
		outdir + "/stats/reads.bxstats.html",
		outdir + "/stats/samtools_stats/bwa.stats.html",
		outdir + "/stats/samtools_flagstat/bwa.flagstat.html"
	message:
		"Read mapping completed!"
	default_target: True

rule link_genome:
	input:
		genomefile
	output: 
		f"Assembly/{bn}"
	message: 
		"Symlinking {input}"
	shell: 
		"ln -sr {input} {output}"


rule faidx_genome:
    input: 
        f"Assembly/{bn}"
    output: 
        f"Assembly/{bn}.fai"
    message:
        "Indexing {input}"
    log:
        f"Assembly/{bn}.faidx.log"
    shell: 
        """
        samtools faidx --fai-idx {output} {input} 2> {log}
        """

rule index_bwa_genome:
    input: 
        f"Assembly/{bn}"
    output: 
        multiext(f"Assembly/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    message:
        "Indexing {input}"
    log:
        f"Assembly/{bn}.idx.log"
    shell: 
        """
        bwa index {input} 2> {log}
        """

rule make_genome_windows:
	input:
		f"Assembly/{bn}.fai"
	output: 
		f"Assembly/{bn}.bed"
	message: 
		"Creating BED intervals from {input}"
	shell: 
		"""
		makewindows.py -i {input} -w 10000 -o {output}
		"""

rule count_beadtags:
	input:
		forward_reads = seq_dir + "/{sample}" + f".{Rsep[0]}.{fqext}",
		reverse_reads = seq_dir + "/{sample}" + f".{Rsep[1]}.{fqext}"
	output: 
		counts = temp(outdir + "/stats/bxcount/{sample}.ema-ncnt"),
		logs   = temp(outdir + "/stats/bxcount/{sample}.count.log")
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message:
		"Counting barcode frequency: {wildcards.sample}"
	params:
		prefix = lambda wc: outdir + "/stats/bxcount/" + wc.get("sample")
	threads: 1
	shell:
		"seqfu interleave -1 {input.forward_reads} -2 {input.reverse_reads} | ema-h count -p -o {params} 2> {output.logs}"

rule beadtag_summary:
	input: 
		countlog = expand(outdir + "/stats/bxcount/{sample}.count.log", sample = samplenames)
	output:
		outdir + "/stats/reads.bxstats.html"
	message:
		"Creating sample barcode validation report"
	script:
		"reportBxCount.Rmd"

rule align:
	input:
		forward_reads = seq_dir + "/{sample}" + f".{Rsep[0]}.{fqext}",
		reverse_reads = seq_dir + "/{sample}" + f".{Rsep[1]}.{fqext}",
		genome 		  = f"Assembly/{bn}",
		genome_idx 	  = multiext(f"Assembly/{bn}", ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
	output:  
		bam    = temp(outdir + "/{sample}.sort.bam"),
		tmpdir = temp(directory(outdir + "/{sample}"))
	log:
		outdir + "/logs/{sample}.log"
	message:
		"Aligning sequences: {wildcards.sample}"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	benchmark:
		"Benchmark/Mapping/bwa/align.{sample}.txt"
	params: 
		quality = config["quality"],
		extra   = extra
	threads:
		8
	shell:
		"""
		mkdir -p Align/bwa/{wildcards.sample}
		BWA_THREADS=$(( {threads} - 2 ))
		bwa mem -C -t $BWA_THREADS {params.extra} -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" {input.genome} {input.forward_reads} {input.reverse_reads} 2> {log} |
		samtools view -h -q {params.quality} | 
		samtools sort -T Align/bwa/{wildcards.sample} --reference {input.genome} -O bam -l 0 -m 4G -o {output.bam} 2> /dev/null
		"""

rule mark_duplicates:
	input:
		outdir + "/{sample}.sort.bam"
	output:
		bam = outdir + "/{sample}.bam",
		bai = outdir + "/{sample}.bam.bai"
	log:
		outdir + "/logs/{sample}.markdup.log"
	message:
		f"Marking duplicates: " + "{wildcards.sample}"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	benchmark:
		"Benchmark/Mapping/bwa/markdup.{sample}.txt"
	threads: 
		4
	shell:
		"sambamba markdup -t {threads} -l 0 {input} {output.bam} 2> {log}"

rule alignment_coverage:
	input: 
		bed = f"Assembly/{bn}.bed",
		bam = outdir + "/{sample}.bam"
	output: 
		outdir + "/stats/coverage/data/{sample}.cov.gz"
	message:
		"Calculating genomic coverage: {wildcards.sample}"
	threads: 
		2
	shell:
		"samtools bedcov -c {input} | gzip > {output}"

rule coverage_report:
	input:
		outdir + "/stats/coverage/data/{sample}.cov.gz"
	output:
		outdir + "/stats/coverage/{sample}.cov.html"
	message:
		"Summarizing alignment coverage: {wildcards.sample}"
	script:
		"reportBwaGencov.Rmd"

rule BEDconvert:
	input:
		bam = outdir + "/{sample}.bam"
	output: 
		unfilt = temp(outdir + "/bedfiles/{sample}.bed"),
		bx     = temp(outdir + "/bedfiles/{sample}.bx.bed")
	message:
		"Converting to BED format: {wildcards.sample}"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	params: lambda wc: outdir + "/align/" + wc.get("sample") + "/" + wc.get("sample") + ".bed"
	threads: 1
	shell:
		"""
		writeBED.pl {input} {output.unfilt}
		awk '!($4~/A00|B00|C00|D00/)' {output.unfilt} > {output.bx}
		"""

rule BX_stats:
	input:
		bedfile  = outdir + "/bedfiles/{sample}.bx.bed"
	output:	
		molsize  = outdir + "/stats/moleculesize/{sample}.molsize.gz",
		molhist  = outdir + "/stats/moleculesize/{sample}.molsize.hist",
		readsper = outdir + "/stats/readsperbx/{sample}.readsperbx"
	message: 
		"Calculating molecule size, reads per molecule: {wildcards.sample}"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	threads: 1
	shell:
		"""
		cut -f10 {input} | datamash -s groupby 1 count 1 | sort -k 1 -n > {output.readsper}
		awk '{{ print $1"\\t"$2"\\t"$3"\\t"$3-$2"\\t"$4"\\t"$10 }}' {input} | sort -k 4 -n | gzip > {output.molsize}
		zcat {output.molsize} | cut -f4 | datamash bin:1000 1 | datamash -s groupby 1 count 1 | sort -k 1 -n > {output.molhist}
		"""

rule alignment_stats:
	input:
		bam      = outdir + "/{sample}.bam",
		bai      = outdir + "/{sample}.bam.bai"
	output: 
		stats    = outdir + "/stats/samtools_stats/{sample}.stats",
		flagstat = outdir + "/stats/samtools_flagstat/{sample}.flagstat"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message:
		"Calculating alignment stats: {wildcards.sample}"
	benchmark:
		"Benchmark/Mapping/bwa/stats.{sample}.txt"
	shell:
		"""
		samtools stats {input.bam} > {output.stats}
		samtools flagstat {input.bam} > {output.flagstat}
		"""

rule samtools_reports:
	input: 
		expand(outdir + "/stats/samtools_{ext}/{sample}.{ext}", sample = samplenames, ext = ["stats", "flagstat"])
	output: 
		stats    = outdir + "/stats/samtools_stats/bwa.stats.html",
		flagstat = outdir + "/stats/samtools_flagstat/bwa.flagstat.html"
	message:
		"Summarizing samtools stats and flagstats"
	shell:
		"""
		multiqc Align/bwa/stats/samtools_stats    --force --quiet --no-data-dir --filename {output.stats} 2> /dev/null
		multiqc Align/bwa/stats/samtools_flagstat --force --quiet --no-data-dir --filename {output.flagstat} 2> /dev/null
		"""