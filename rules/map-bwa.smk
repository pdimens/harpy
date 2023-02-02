# user specified configs
seq_dir = config["seq_directory"]
genomefile = config["genome_file"]
# Received from the harpy wrapper
Rsep = config["Rsep"]
fqext = config["fqext"]
samplenames = config["samplenames"]
extra = config["extra"]

rule create_reports:
	input: 
		expand("ReadMapping/align/{sample}.bam", sample = samplenames),
		expand("ReadMapping/align/{ext}/{sample}.{ext}", sample = samplenames, ext = ["stats", "flagstat"])
	output: 
		stats = "ReadMapping/alignment.stats.html",
		flagstat = "ReadMapping/alignment.flagstat.html"
	message: "Read mapping completed!\nAlignment reports:\n{output.stats}\n{output.flagstat}"
	benchmark: "Benchmark/Mapping/bwa/reports.txt"
	default_target: True
	shell:
		"""
		multiqc ReadMapping/align/stats --force --quiet --no-data-dir --filename {output.stats} 2> /dev/null
		multiqc ReadMapping/align/flagstat --force --quiet --no-data-dir --filename {output.flagstat} 2> /dev/null
		"""

rule index_genome:
	input: genomefile
	output: multiext(genomefile, ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
	message: "Indexing {input} prior to read mapping"
	benchmark: "Benchmark/Mapping/bwa/genoindex.txt"
	shell: 
		"""
		bwa index {input}
		samtools faidx {input}
		"""

rule align_bwa:
	input:
		forward_reads = seq_dir + "/{sample}" + Rsep + "1." + fqext,
		reverse_reads = seq_dir + "/{sample}" + Rsep + "2." + fqext,
		genome = genomefile,
		genome_idx = multiext(genomefile, ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
	output:  pipe("ReadMapping/align/{sample}.sam")
	log: "ReadMapping/align/logs/{sample}.log"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message: "Mapping {wildcards.sample} reads onto {input.genome} using BWA"
	benchmark: "Benchmark/Mapping/bwa/align.{sample}.txt"
	params: 
		extra = extra
	log: "ReadMapping/count/logs/{sample}.count.log"
	threads: 3
	shell:
		"""
		bwa mem -p -C -t {threads} {params} -M -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" {input.genome} {input.forward_reads} {input.reverse_reads} > {output} 2> {log}
		"""

rule sort_alignments:
	input: "ReadMapping/align/{sample}.sam"
	output: 
		bam = temp("ReadMapping/align/{sample}.sort.bam")
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message: "Sorting {wildcards.sample} alignments"
	benchmark: "Benchmark/Mapping/bwa/sort.{sample}.txt"
	threads: 1
	shell:
		"""
		samtools sort -@ {threads} -O bam -l 0 -m 4G -o {output} {input}
		"""

rule mark_duplicates:
	input: "ReadMapping/align/{sample}.sort.bam"
	output: 
		bam = "ReadMapping/align/{sample}.bam",
		bai = "ReadMapping/align/{sample}.bam.bai"
	log: "ReadMapping/align/log/{sample}.markdup.nobarcode.log"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message: "Marking duplicates: {wildcards.sample}"
	benchmark: "Benchmark/Mapping/bwa/markdup.{sample}.txt"
	threads: 4
	shell:
		"""
		sambamba markdup -t {threads} -l 0 {input} {output.bam} 2> {log}
		samtools index {output.bam}
		"""

rule alignment_stats:
	input:
		bam = "ReadMapping/align/{sample}.bam",
		bai = "ReadMapping/align/{sample}.bam.bai"
	output: 
		stats = "ReadMapping/align/stats/{sample}.stats",
		flagstat = "ReadMapping/align/flagstat/{sample}.flagstat"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message: "Calculating alignment stats: {wildcards.sample}"
	benchmark: "Benchmark/Mapping/bwa/stats.{sample}.txt"
	threads: 1
		"""
		samtools stats {input.bam} > {output.stats}
		samtools flagstat {input.bam} > {output.flagstat}
		"""