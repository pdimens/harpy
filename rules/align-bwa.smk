import os

seq_dir = config["seq_directory"]
genomefile = config["genomefile"]
Rsep = config["Rsep"]
fqext = config["fqext"]
samplenames = config["samplenames"]
extra = config.get("extra", "") 

bn = os.path.basename(genomefile)

rule create_reports:
	input: 
		expand("Alignments/bwa/{sample}.bam", sample = samplenames),
		expand("Alignments/bwa/stats/coverage/{sample}.cov.html", sample = samplenames),
		"Alignments/bwa/stats/samtools_stats/bwa.stats.html",
		"Alignments/bwa/stats/samtools_flagstat/bwa.flagstat.html"
	message: "Read mapping completed!"
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


rule index_genome:
	input:
		f"Assembly/{bn}"
	output: 
		multiext(f"Assembly/{bn}", ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb", ".bed")
	message:
		"Indexing {input}"
	log: 
		f"Assembly/{bn}.idx.log"
	shell: 
		"""
		bwa index {input} 2> {log}
		samtools faidx --fai-idx {input}.fai {input} 2>> {log}
		utilities/makewindows.py -i {input}.fai -w 10000 -o {input}.bed
		"""

rule align:
	input:
		forward_reads = seq_dir + "/{sample}" + f".{Rsep[0]}.{fqext}",
		reverse_reads = seq_dir + "/{sample}" + f".{Rsep[1]}.{fqext}",
		genome = f"Assembly/{bn}",
		genome_idx = multiext(f"Assembly/{bn}", ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
	output:  
		bam = temp("Alignments/bwa/{sample}.sort.bam"),
		tmpdir = temp(directory("Alignments/bwa/{sample}"))
	log: "Alignments/bwa/logs/{sample}.log"
	message: "Aligning sequences: {wildcards.sample}"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	benchmark: "Benchmark/Mapping/bwa/align.{sample}.txt"
	params: 
		quality = config["quality"],
		extra = extra
	threads: 8
	shell:
		"""
		mkdir -p Alignments/bwa/{wildcards.sample}
		BWA_THREADS=$(( {threads} - 2 ))
		bwa mem -C -t $BWA_THREADS {params.extra} -M -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" {input.genome} {input.forward_reads} {input.reverse_reads} 2> {log} |
		samtools view -h -q {params.quality} | 
		samtools sort -T Alignments/bwa/{wildcards.sample} --reference {input.genome} -O bam -l 0 -m 4G -o {output.bam} 2> /dev/null
		"""

#rule sort_alignments:
#	input: 
#		sam = "Alignments/bwa/{sample}.sam",
#		asm = f"Assembly/{genomefile}"
#	output: temp("Alignments/bwa/{sample}.sort.bam")
#	wildcard_constraints:
#		sample = "[a-zA-Z0-9_-]*"
#	message: "Sorting {wildcards.sample} alignments"
#	benchmark: "Benchmark/Mapping/bwa/sort.{sample}.txt"
#	threads: 2
#	shell:
#		"""
#		samtools sort --threads {threads} --reference {input.asm} -O bam -m 4G -o {output} {input.sam}
#		"""

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
	threads: 4
	shell:
		"sambamba markdup -t {threads} -l 0 {input} {output.bam} 2> {log}"

rule alignment_coverage:
	input: 
		bed = f"Assembly/{bn}.bed",
		bam = "Alignments/bwa/{sample}.bam"
	output: "Alignments/bwa/stats/coverage/data/{sample}.cov.gz"
	message: "Calculating genomic coverage: {wildcards.sample}"
	threads: 2
	shell:
		"samtools bedcov -c {input} | gzip > {output}"

rule coverage_report:
	input:
		gencov = "Alignments/bwa/stats/coverage/data/{sample}.cov.gz"
	output:
		"Alignments/bwa/stats/coverage/{sample}.cov.html"
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