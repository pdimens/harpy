# user specified configs
seq_dir = config["seq_directory"]
nbins = config["EMA_bins"]
genomefile = config["genomefile"]
Rsep = config["Rsep"]
fqext = config["fqext"]
samplenames = config["samplenames"]
extra = config.get("extra", "") 

rule create_reports:
	input: 
		expand("ReadMapping/align/{sample}.bam", sample = samplenames),
		expand("ReadMapping/align/{sample}.bam.bai", sample = samplenames),
		expand("ReadMapping/align/{ext}/{sample}.{ext}", sample = samplenames, ext = ["stats", "flagstat"]),
		"ReadMapping/count/Beadtag.report"
	output: 
		stats = "ReadMapping/alignment.stats.html",
		flagstat = "ReadMapping/alignment.flagstat.html"
	message: "Read mapping completed!\nAlignment reports:\n{output.stats}\n{output.flagstat}"
	benchmark: "Benchmark/Mapping/ema/report.txt"
	default_target: True
	shell:
		"""
		multiqc ReadMapping/align/stats --force --quiet --no-data-dir --filename {output.stats} 2> /dev/null
		multiqc ReadMapping/align/flagstat --force --quiet --no-data-dir --filename {output.flagstat} 2> /dev/null
		"""

rule index_genome:
	input: genomefile
	output: multiext(genomefile, ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
	message: "Indexing {input}"
	benchmark: "Benchmark/Mapping/genoindex.txt"
	shell: 
		"""
		bwa index {input}
		samtools faidx {input}
		"""

rule count_beadtags:
	input:
		forward_reads = seq_dir + "/{sample}" + f".{Rsep[0]}.{fqext}",
		reverse_reads = seq_dir + "/{sample}" + f".{Rsep[1]}.{fqext}"
	output: 
		counts = "ReadMapping/count/{sample}.ema-ncnt",
		logs = temp("ReadMapping/count/logs/{sample}.count.log")
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message: "Counting barcode frequency: {wildcards.sample}"
	benchmark: "Benchmark/Mapping/ema/Count.{sample}.txt"
	params:
		prefix = lambda wc: "ReadMapping/count/" + wc.get("sample")
	threads: 1
	shell:
		"""
		emaInterleave {input.forward_reads} {input.reverse_reads} | ema-h count -p -o {params} 2> {output.logs}
		"""

rule beadtag_summary:
	input: expand("ReadMapping/count/logs/{sample}.count.log", sample = samplenames)
	output: "ReadMapping/count/Beadtag.report"
	message: "Creating sample barcode validation report"
	benchmark: "Benchmark/Mapping/ema/beadtagsummary.txt"
	run:
		import os
		with open(output[0], "w") as outfile:
			outfile.write('{0!s} {1!s} {2!s} {3!s}\n'.format("Sample", "BarcodeOK", "BarcodeTotal", "PercentTotal"))
			for i in input:
				with open(i) as f:
					samplename = os.path.basename(i).split(".count.log")[0]
					bc_line = [line for line in f.read().splitlines() if 'OK barcode:' in line][0].split(" ")
					bc_len = len(bc_line)
					bc_ok = bc_line[bc_len - 4]
					bc_total = bc_line[bc_len - 1]
					bc_percent = int(bc_ok.replace(',', '')) / int(bc_total.replace(',', '')) * 100
					outfile.write('{0!s} {1!s} {2!s} {3!s}\n'.format(samplename, bc_ok, bc_total, round(bc_percent, 2)))

rule preprocess_ema:
	input: 
		forward_reads = seq_dir + "/{sample}" + f".{Rsep[0]}.{fqext}",
		reverse_reads = seq_dir + "/{sample}" + f".{Rsep[1]}.{fqext}",
		emacounts = "ReadMapping/count/{sample}.ema-ncnt"
	output: 
		bins = temp(expand("ReadMapping/preproc/{{sample}}/ema-bin-{bin}", bin = ["%03d" % i for i in range(nbins)])),
		unbarcoded = temp("ReadMapping/preproc/{sample}/ema-nobc")
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	log: "ReadMapping/preproc/logs/{sample}.preproc.log"
	message: "Preprocessing for EMA mapping: {wildcards.sample}"
	benchmark: "Benchmark/Mapping/ema/Preproc.{sample}.txt"
	threads: 2
	params:
		outdir = lambda wc: "ReadMapping/preproc/" + wc.get("sample"),
		bins = nbins
	shell:
		"""
		emaInterleave {input.forward_reads} {input.reverse_reads} | ema-h preproc -p -n {params.bins} -t {threads} -o {params.outdir} {input.emacounts} 2>&1 | cat - > {log}
		"""

rule align_ema:
	input:
		readbin = "ReadMapping/preproc/{sample}/ema-bin-{bin}",
		genome = genomefile,
		genome_idx = multiext(genomefile, ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
	output: pipe("ReadMapping/align/{sample}/{sample}.{bin}.sam")
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
		reads = "ReadMapping/preproc/{sample}/ema-nobc",
		genome = genomefile,
		genome_idx = multiext(genomefile, ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
	output: 
		samfile = pipe("ReadMapping/align/{sample}/{sample}.nobarcode.sam")
	benchmark: "Benchmark/Mapping/ema/bwaAlign.{sample}.txt"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message: "Mapping unbarcoded reads onto {input.genome}: {wildcards.sample}"
	threads: 2
	shell:
		"""
		bwa mem -p -t {threads} -M -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" {input.genome} {input.reads} > {output} 2> /dev/null
		"""

rule sort_ema:
	input: "ReadMapping/align/{sample}/{sample}.{emabin}.sam"
	output: temp("ReadMapping/align/{sample}/{sample}.{emabin}.bam")
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*",
		emabin = "[0-9]*"
	message: "Sorting alignments: {wildcards.sample}-{wildcards.emabin}"
	benchmark: "Benchmark/Mapping/ema/Sort.{sample}.{emabin}.txt"
	threads: 1
	shell: 
		"""
		samtools sort -@ {threads} -O bam -l 0 -m 4G -o {output} {input}
		"""

rule sort_nobarcode:
	input: "ReadMapping/align/{sample}/{sample}.nobarcode.sam"
	output: temp("ReadMapping/align/{sample}/{sample}.nobarcode.bam.tmp")
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message: "Sorting unbarcoded alignments: {wildcards.sample}"
	benchmark: "Benchmark/Mapping/ema/bwaSort.{sample}.txt"
	threads: 2
	shell:
		"""
		samtools sort -@ {threads} -O bam -l 0 -m 4G -o {output} {input}
		"""    

rule markduplicates:
	input: "ReadMapping/align/{sample}/{sample}.nobarcode.bam.tmp"
	output: 
		bam = temp("ReadMapping/align/{sample}/{sample}.nobarcode.bam"),
		bai = temp("ReadMapping/align/{sample}/{sample}.nobarcode.bam.bai")
	log: 
		mdlog = "ReadMapping/align/log/{sample}.markdup.nobarcode.log",
		stats = "ReadMapping/align/stats/{sample}.nobarcode.stats",
		flagstat = "ReadMapping/align/flagstat/{sample}.nobarcode.flagstat"
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
		aln_barcoded = expand("ReadMapping/align/{{sample}}/{{sample}}.{bin}.bam", bin = ["%03d" % i for i in range(nbins)]),
	output: temp("ReadMapping/align/{sample}/{sample}.barcoded.nosort.bam"),
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message: "Merging barcoded alignments: {wildcards.sample}"
	benchmark: "Benchmark/Mapping/ema/merge.{sample}.txt"
	threads: 10
	shell:
		"""
		sambamba merge -t {threads} {output} {input}
		"""	

rule sort_barcoded:
	input:
		bam = "ReadMapping/align/{sample}/{sample}.barcoded.nosort.bam",
		genome = genomefile
	output: "ReadMapping/align/{sample}/{sample}.barcoded.bam"
	message: "Sorting merged barcoded alignments: {wildcards.sample}"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	threads: 2
	shell:
		"""
		samtools sort -@ {threads} -O bam --reference {input.genome} -l 0 -m 4G -o {output} {input.bam}
		"""

rule index_mergedbarcoded:
	input: "ReadMapping/align/{sample}/{sample}.barcoded.bam"
	output: temp("ReadMapping/align/{sample}/{sample}.barcoded.bam.bai")
	log:
		stats = "ReadMapping/align/stats/{sample}.barcoded.stats",
		flagstat = "ReadMapping/align/flagstat/{sample}.barcoded.flagstat"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message: "Indexing merged barcoded alignemnts: {wildcards.sample}"
	benchmark: "Benchmark/Mapping/ema/indexmerge.{sample}.txt"
	shell:
		"""
		sambamba index {input} {output}
		samtools stats {input} > {log.stats}
		samtools flagstat {input} > {log.flagstat}
		"""

rule merge_alignments:
	input:
		aln_barcoded = "ReadMapping/align/{sample}/{sample}.barcoded.bam",
		aln_nobarcode = "ReadMapping/align/{sample}/{sample}.nobarcode.bam",
		idx_barcoded = "ReadMapping/align/{sample}/{sample}.barcoded.bam.bai",
		idx_nobarcode = "ReadMapping/align/{sample}/{sample}.nobarcode.bam.bai"
	output: 
		bam = "ReadMapping/align/{sample}.bam"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message: "Merging all alignments: {wildcards.sample}"
	benchmark: "Benchmark/Mapping/ema/mergebc_nobc.{sample}.txt"
	threads: 10
	shell:
		"""
		sambamba merge -t {threads} {output.bam} {input.aln_barcoded} {input.aln_nobarcode}
		"""

rule index_alignments:
	input: "ReadMapping/align/{sample}.bam"
	output: "ReadMapping/align/{sample}.bam.bai"
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
		bam = "ReadMapping/align/{sample}.bam",
		bai = "ReadMapping/align/{sample}.bam.bai"
	output:
		stats = report("ReadMapping/align/stats/{sample}.stats", category="{sample}", subcategory="All Aligments", labels={"Metric": "stats"}),
		flagstat = report("ReadMapping/align/flagstat/{sample}.flagstat", category="{sample}", subcategory="All Aligments", labels={"Metric": "flagstat"})
	message: "Calculating alignment stats: {wildcards.sample}"
	benchmark: "Benchmark/Mapping/ema/Mergedstats.{sample}.txt"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	shell:
		"""
		samtools stats {input.bam} > {output.stats}
		samtools flagstat {input.bam} > {output.flagstat}
		"""