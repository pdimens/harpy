import os

seq_dir = config["seq_directory"]
genomefile = config["genomefile"]
Rsep = config["Rsep"]
fqext = config["fqext"]
samplenames = config["samplenames"]
extra = config.get("extra", "") 

bn = os.path.basename(genomefile)
shell("mkdir -p Assembly")
if not os.path.exists(f"Assembly/{bn}"):
    shell(f"ln -sr {genomefile} Assembly/{bn}")

rule create_reports:
	input: 
		expand("ReadMapping/bwa/{sample}.bam", sample = samplenames),
		expand("ReadMapping/bwa/stats/samtools_{ext}/{sample}.{ext}", sample = samplenames, ext = ["stats", "flagstat"]),
		expand("ReadMapping/bwa/stats/coverage/{sample}.gencov", sample = samplenames)
	output: 
		stats =    "ReadMapping/bwa/stats/samtools_stats/bwa.stats.html",
		flagstat = "ReadMapping/bwa/stats/samtools_flagstat/bwa.flagstat.html"
	message: "Read mapping completed!\nAlignment reports:\n{output.stats}\n{output.flagstat}"
	benchmark: "Benchmark/Mapping/bwa/reports.txt"
	default_target: True
	shell:
		"""
		multiqc ReadMapping/bwa/stats/samtools_stats    --force --quiet --no-data-dir --filename {output.stats} 2> /dev/null
		multiqc ReadMapping/bwa/stats/samtools_flagstat --force --quiet --no-data-dir --filename {output.flagstat} 2> /dev/null
		"""

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
	input: 
		sam = "ReadMapping/bwa/{sample}.sam",
		asm = f"Assembly/{genomefile}"
	output: temp("ReadMapping/bwa/{sample}.sort.bam")
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
	input: "ReadMapping/bwa/{sample}.sort.bam"
	output:
		bam = "ReadMapping/bwa/{sample}.bam",
		bai = "ReadMapping/bwa/{sample}.bam.bai"
	log: "ReadMapping/bwa/logs/{sample}.markdup.log"
	message: f"Marking duplicates: " + "{wildcards.sample}"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	benchmark: "Benchmark/Mapping/bwa/markdup.{sample}.txt"
	params:
		rootname = "ReadMapping/bwa/{sample}"
	threads: 4
	shell:
		"""
		sambamba markdup -t {threads} -l 0 {input} {output.bam} 2> {log}
		"""

#	if [[ "{params.bx}" == "False" ]]; then
#	else
#		samtools collate -o {params.rootname}.collate.bam {input} 2> {log}
#		samtools fixmate -m {params.rootname}.collate.bam {params.rootname}.fixmate.bam 2>> {log}
#		rm {params.rootname}.collate.bam
#		samtools sort -O bam {params.rootname}.fixmate.bam > {params.rootname}.fixsort.bam 2>> {log}
#		rm {params.rootname}.fixmate.bam
#		samtools markdup --barcode-tag BX {params.rootname}.fixsort.bam {output.bam} 2>> {log}
#		rm {params.rootname}.fixsort.bam
#		samtools index {output.bam} 2> {log}
#	fi

rule genome_coords:
	input: f"Assembly/{genomefile}.fai"
	output: f"Assembly/{genomefile}.bed"
	message: "Creating BED file of genomic coordinates"
	threads: 1
	shell:
		"""
		awk 'BEGIN{{FS=OFS="\\t"}} {{print $1,$2}}' {input} > {output}
		"""

rule BEDconvert:
	input: "ReadMapping/bwa/{sample}.bam"
	output: temp("ReadMapping/bwa/bedfiles/{sample}.bed")
	message: "Converting to BED format: {wildcards.sample}"
	shell:
		"bedtools bamtobed -i {input} > {output}"

rule genome_coverage:
	input:
		geno = f"Assembly/{genomefile}.bed",
		bed = "ReadMapping/bwa/bedfiles/{sample}.bed"
	output: 
		"ReadMapping/bwa/stats/coverage/{sample}.gencov"
	message: 
		"Calculating genomic coverage of alignments: {wildcards.sample}"
	shell:
		"""
		bedtools genomecov -i {input.bed} -g {input.geno} > {output}
		"""

rule alignment_stats:
	input:
		bam = "ReadMapping/bwa/{sample}.bam",
		bai = "ReadMapping/bwa/{sample}.bam.bai"
	output: 
		stats = "ReadMapping/bwa/stats/samtools_stats/{sample}.stats",
		flagstat = "ReadMapping/bwa/stats/samtools_flagstat/{sample}.flagstat"
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