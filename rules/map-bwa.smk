import os

# user specified configs
seq_dir = config["seq_directory"]
genomefile = config["genome_file"]

# this identifies whether .fastq.gz or .fq.gz is used as the file extension
fastqlist = [i for i in os.listdir(seq_dir) if i.endswith('.fastq.gz')]
fqext = "fq.gz" if not fastqlist else "fastq.gz"

Rlist = [i for i in os.listdir(seq_dir) if i.endswith('.R1.' + fqext)]
Rsep = "_R" if not Rlist else ".R"
fullext = Rsep + "1." + fqext
samplenames = set([i.split(fullext)[0] for i in os.listdir(seq_dir) if i.endswith(fullext)])
assert(len(samplenames) > 0), "No alignment (.bam) files found in " + bam_dir

rule create_reports:
	input: 
		expand("ReadMapping/align/{sample}.{ext}", sample = samplenames, ext = ["bam", "stats", "flagstat"])
	output: 
		stats = report("ReadMapping/alignment.stats.html", caption = "Samtools stats alignment metrics"),
		flagstat = report("ReadMapping/alignment.flagstat.html", caption = "Samtools flagstat alignment metrics")
	message: "Read mapping completed!\nAlignment reports:\n{output.stats}\n{output.flagstat}"
	default_target: True
	shell:
		"""
		multiqc ReadMapping/align/stats --force --quiet --filename {output.stats}
		multiqc ReadMapping/align/flagstat --force --quiet --filename {output.flagstat}
		"""

rule index_genome:
	input: genomefile
	output: multiext(genomefile, ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
	message: "Indexing {input} prior to read mapping"
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
	log: "ReadMapping/count/logs/{sample}.count.log"
	threads: 3
	shell:
		"""
		bwa mem -p -C -t {threads} -M -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" {input.genome} {input.forward_reads} {input.reverse_reads} > {output} 2> {log}
		"""

rule sort_alignments:
	input: "ReadMapping/align/{sample}.sam"
	output: 
		bam = temp("ReadMapping/align/{sample}.bam.tmp")
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message: "Sorting {wildcards.sample} alignments"
	threads: 1
	shell:
		"""
		samtools sort -@ {threads} -O bam -l 0 -m 4G -o {output} {input}
		"""

rule mark_duplicates:
	input: "ReadMapping/align/{sample}.bam.tmp"
	output: 
		bam = "ReadMapping/align/{sample}.bam",
		bai = "ReadMapping/align/{sample}.bam.bai",
		stats = "ReadMapping/align/stats/{sample}.stats",
		flagstat = "ReadMapping/align/flagstat/{sample}.flagstat"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message: "Marking duplicates with Sambamba for {wildcards.sample} alignments and calculating alignment stats"
	threads: 4
	shell:
		"""
		sambamba markdup -t {threads} -l 0 {input} {output}
		samtools index {output.bam}
		samtools stats {output.bam} > {output.stats}
		samtools flagstat {output.bam} > {output.flagstat}
		"""
 
 