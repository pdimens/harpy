import subprocess

seq_dir = config["seq_directory"]
genomefile = config["genomefile"]
Rsep = config["Rsep"]
fqext = config["fqext"]
samplenames = config["samplenames"]
BXmarkdup = config["BXmarkdup"]
extra = config.get("extra", "") 
txt = " using BX barcodes" if BXmarkdup else ""

rule create_reports:
	input: 
		expand("ReadMapping/bwa/{sample}.bam", sample = samplenames),
		expand("ReadMapping/bwa/{ext}/{sample}.{ext}", sample = samplenames, ext = ["stats", "flagstat"]),
		expand("ReadMapping/bwa/coverage/{sample}.gencov", sample = samplenames)
	output: 
		stats = "ReadMapping/alignment.stats.html",
		flagstat = "ReadMapping/alignment.flagstat.html"
	message: "Read mapping completed!\nAlignment reports:\n{output.stats}\n{output.flagstat}"
	benchmark: "Benchmark/Mapping/bwa/reports.txt"
	default_target: True
	shell:
		"""
		multiqc ReadMapping/bwa/stats --force --quiet --no-data-dir --filename {output.stats} 2> /dev/null
		multiqc ReadMapping/bwa/flagstat --force --quiet --no-data-dir --filename {output.flagstat} 2> /dev/null
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
		forward_reads = seq_dir + "/{sample}" + f".{Rsep[0]}.{fqext}",
		reverse_reads = seq_dir + "/{sample}" + f".{Rsep[1]}.{fqext}",
		genome = genomefile,
		genome_idx = multiext(genomefile, ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
	output:  pipe("ReadMapping/bwa/{sample}.sam")
	log: "ReadMapping/bwa/logs/{sample}.log"
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
	input: "ReadMapping/bwa/{sample}.sam"
	output: temp("ReadMapping/bwa/{sample}.sort.bam")
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message: "Sorting {wildcards.sample} alignments"
	benchmark: "Benchmark/Mapping/bwa/sort.{sample}.txt"
	threads: 1
	shell:
		"""
		samtools sort --threads {threads} -O bam -l 0 -m 4G -o {output} {input}
		"""

rule mark_duplicates:
	input: "ReadMapping/bwa/{sample}.sort.bam"
	output: 
		bam = "ReadMapping/bwa/{sample}.bam",
		bai = "ReadMapping/bwa/{sample}.bam.bai"
	log: "ReadMapping/bwa/log/{sample}.markdup.log"
	message: f"Marking duplicates{txt}: " + "{wildcards.sample}"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	benchmark: "Benchmark/Mapping/bwa/markdup.{sample}.txt"
	threads: 4
	run:
		if BXmarkdup:
			subprocess.run(f"samtools markdup --threads {threads} --barcode-tag BX {input[0]} {output.bam} 2> {log[0]}".split())
		else:
			subprocess.run(f"sambamba markdup -t {threads} -l 0 {input[0]} {output.bam} 2> {log[0]}".split())
#
#rule genome_coords:
#	input: genomefile + ".fai"
#	output: genomefile + ".bed"
#	message: "Creating BED file of genomic coordinates"
#	threads: 1
#	shell:
#		"""
#		awk 'BEGIN{{FS=OFS="\\t"}} {{print $1,$2}}' {input} > {output}
#		"""
#
#rule BEDconvert:
#	input: "ReadMapping/bwa/{sample}.bam"
#	output: temp("ReadMapping/bedfiles/{sample}.bed")
#	message: "Converting to BED format: {wildcards.sample}"
#	shell:
#		"bedtools bamtobed -i {input}"
#
#rule genome_coverage:
#	input:
#		geno = genomefile + ".bed",
#		bed = "ReadMapping/bedfiles/{sample}.bed"
#	output: 
#		"ReadMapping/bwa/coverage/{sample}.gencov"
#	message: 
#		"Calculating genomic coverage of alignments: {wildcards.sample}"
#	shell:
#		"""
#		bedtools genomecov -i {input.bed} -g {input.geno} > {output}
#		"""
#
#rule alignment_stats:
#	input:
#		bam = "ReadMapping/bwa/{sample}.bam",
#		bai = "ReadMapping/bwa/{sample}.bam.bai"
#	output: 
#		stats = "ReadMapping/bwa/stats/{sample}.stats",
#		flagstat = "ReadMapping/bwa/flagstat/{sample}.flagstat"
#	wildcard_constraints:
#		sample = "[a-zA-Z0-9_-]*"
#	message: "Calculating alignment stats: {wildcards.sample}"
#	benchmark: "Benchmark/Mapping/bwa/stats.{sample}.txt"
#	threads: 1
#		"""
#		samtools stats {input.bam} > {output.stats}
#		samtools flagstat {input.bam} > {output.flagstat}
#		"""