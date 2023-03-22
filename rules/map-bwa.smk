import os

seq_dir = config["seq_directory"]
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
		expand("Alignments/bwa/{sample}.bam", sample = samplenames),
		expand("Alignments/bwa/stats/coverage/{sample}.gencov.html", sample = samplenames),
		"Alignments/bwa/stats/samtools_stats/bwa.stats.html",
		"Alignments/bwa/stats/samtools_flagstat/bwa.flagstat.html"
	message: "Read mapping completed!"
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

rule align_bwa:
	input:
		forward_reads = seq_dir + "/{sample}" + f".{Rsep[0]}.{fqext}",
		reverse_reads = seq_dir + "/{sample}" + f".{Rsep[1]}.{fqext}",
		genome = f"Assembly/{genomefile}",
		genome_idx = multiext(f"Assembly/{genomefile}", ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
	output:  pipe("Alignments/bwa/{sample}.sam")
	log: "Alignments/bwa/logs/{sample}.log"
	message: "Mapping onto {input.genome}: {wildcards.sample}"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	benchmark: "Benchmark/Mapping/bwa/align.{sample}.txt"
	params: 
		extra = extra
	threads: 3
	shell:
		"""
		bwa mem -C -t {threads} {params} -M -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" {input.genome} {input.forward_reads} {input.reverse_reads} > {output} 2> {log}
		"""

rule sort_alignments:
	input: 
		sam = "Alignments/bwa/{sample}.sam",
		asm = f"Assembly/{genomefile}"
	output: temp("Alignments/bwa/{sample}.sort.bam")
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message: "Sorting {wildcards.sample} alignments"
	benchmark: "Benchmark/Mapping/bwa/sort.{sample}.txt"
	threads: 1
	shell:
		"""
		samtools sort --threads {threads} --reference {input.asm} -O bam -l 0 -m 4G -o {output} {input.sam}
		"""

rule mark_duplicates:
	input: "Alignments/bwa/{sample}.sort.bam"
	output:
		bam = "Alignments/bwa/{sample}.bam",
		bai = "Alignments/bwa/{sample}.bam.bai"
	log: "Alignments/bwa/logs/{sample}.markdup.log"
	message: f"Marking duplicates: " + "{wildcards.sample}"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	benchmark: "Benchmark/Mapping/bwa/markdup.{sample}.txt"
	params:
		rootname = "Alignments/bwa/{sample}"
	threads: 4
	shell:
		"""
		sambamba markdup -t {threads} -l 0 {input} {output.bam} 2> {log}
		"""

rule genome_coverage:
	input: "Alignments/bwa/{sample}.bam"
	output: "Alignments/bwa/stats/coverage/data/{sample}.gencov"
	message: "Calculating genomic coverage: {wildcards.sample}"
	threads: 2
	params: bed_prox
	shell:
		"""
		bedtools genomecov -ibam {input} -bg | bedtools merge -c 4 -o sum -d {params} > {output}
		"""

rule gencov_report:
	input:
		gencov = "Alignments/bwa/stats/coverage/data/{sample}.gencov",
		faidx = f"Assembly/{genomefile}.fai"
	output:
		"Alignments/bwa/stats/coverage/{sample}.gencov.html"
	message:
		"Summarizing alignment coverage: {wildcards.sample}"
	script:
		"../utilities/reportGencov.Rmd"

rule alignment_stats:
	input:
		bam = "Alignments/bwa/{sample}.bam",
		bai = "Alignments/bwa/{sample}.bam.bai"
	output: 
		stats = "Alignments/bwa/stats/samtools_stats/{sample}.stats",
		flagstat = "Alignments/bwa/stats/samtools_flagstat/{sample}.flagstat"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message: "Calculating alignment stats: {wildcards.sample}"
	benchmark: "Benchmark/Mapping/bwa/stats.{sample}.txt"
	threads: 1
	shell:
		"""
		samtools stats {input.bam} > {output.stats}
		samtools flagstat {input.bam} > {output.flagstat}
		"""

rule samtools_reports:
	input: 
		expand("Alignments/bwa/stats/samtools_{ext}/{sample}.{ext}", sample = samplenames, ext = ["stats", "flagstat"])
	output: 
		stats =    "Alignments/bwa/stats/samtools_stats/bwa.stats.html",
		flagstat = "Alignments/bwa/stats/samtools_flagstat/bwa.flagstat.html"
	message: "Summarizing samtools stats and flagstats"
	shell:
		"""
		multiqc Alignments/bwa/stats/samtools_stats    --force --quiet --no-data-dir --filename {output.stats} 2> /dev/null
		multiqc Alignments/bwa/stats/samtools_flagstat --force --quiet --no-data-dir --filename {output.flagstat} 2> /dev/null
		"""