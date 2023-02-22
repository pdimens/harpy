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
		expand("ReadMapping/ema/{sample}.bam", sample = samplenames),
		expand("ReadMapping/ema/{sample}.bam.bai", sample = samplenames),
		expand("ReadMapping/ema/{ext}/{sample}.{ext}", sample = samplenames, ext = ["stats", "flagstat"]),
		expand("ReadMapping/bxstats/moleculesize/{sample}.molsize", sample = samplenames),
		expand("ReadMapping/bxstats/moleculesize/{sample}.molsize.hist", sample = samplenames),
		expand("ReadMapping/bxstats/readsperbx/{sample}.readsperbx", sample = samplenames),
		"ReadMapping/ema/count/Beadtag.report"
	output: 
		stats = "ReadMapping/alignment.stats.html",
		flagstat = "ReadMapping/alignment.flagstat.html"
	message: "Read mapping completed!\nAlignment reports:\n\t{output.stats}\n\t{output.flagstat}"
	benchmark: "Benchmark/Mapping/ema/report.txt"
	default_target: True
	shell:
		"""
		multiqc ReadMapping/ema/stats --force --quiet --no-data-dir --filename {output.stats} 2> /dev/null
		multiqc ReadMapping/ema/flagstat --force --quiet --no-data-dir --filename {output.flagstat} 2> /dev/null
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
		counts = "ReadMapping/ema/count/{sample}.ema-ncnt",
		logs = temp("ReadMapping/ema/count/logs/{sample}.count.log")
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message: "Counting barcode frequency: {wildcards.sample}"
	benchmark: "Benchmark/Mapping/ema/Count.{sample}.txt"
	params:
		prefix = lambda wc: "ReadMapping/ema/count/" + wc.get("sample")
	threads: 1
	shell:
		"""
		emaInterleave {input.forward_reads} {input.reverse_reads} | ema-h count -p -o {params} 2> {output.logs}
		"""

rule beadtag_summary:
	input: expand("ReadMapping/ema/count/logs/{sample}.count.log", sample = samplenames)
	output: "ReadMapping/ema/count/Beadtag.report"
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
		emacounts = "ReadMapping/ema/count/{sample}.ema-ncnt"
	output: 
		bins = temp(expand("ReadMapping/ema/preproc/{{sample}}/ema-bin-{bin}", bin = ["%03d" % i for i in range(nbins)])),
		unbarcoded = temp("ReadMapping/ema/preproc/{sample}/ema-nobc")
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	log: "ReadMapping/ema/preproc/logs/{sample}.preproc.log"
	message: "Preprocessing for EMA mapping: {wildcards.sample}"
	benchmark: "Benchmark/Mapping/ema/Preproc.{sample}.txt"
	threads: 2
	params:
		outdir = lambda wc: "ReadMapping/ema/preproc/" + wc.get("sample"),
		bins = nbins
	shell:
		"""
		emaInterleave {input.forward_reads} {input.reverse_reads} | ema-h preproc -p -n {params.bins} -t {threads} -o {params.outdir} {input.emacounts} 2>&1 | cat - > {log}
		"""

rule align_ema:
	input:
		readbin = "ReadMapping/ema/preproc/{sample}/ema-bin-{bin}",
		genome = genomefile,
		genome_idx = multiext(genomefile, ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
	output: pipe("ReadMapping/ema/{sample}/{sample}.{bin}.sam")
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
		reads = "ReadMapping/ema/preproc/{sample}/ema-nobc",
		genome = genomefile,
		genome_idx = multiext(genomefile, ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
	output: 
		samfile = pipe("ReadMapping/ema/{sample}/{sample}.nobarcode.sam")
	benchmark: "Benchmark/Mapping/ema/bwaAlign.{sample}.txt"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message: "Mapping unbarcoded reads onto {input.genome}: {wildcards.sample}"
	threads: 2
	shell:
		"""
		bwa mem -p -t {threads} -C -M -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" {input.genome} {input.reads} > {output} 2> /dev/null
		"""

rule sort_ema:
	input: "ReadMapping/ema/{sample}/{sample}.{emabin}.sam"
	output: temp("ReadMapping/ema/{sample}/{sample}.{emabin}.bam")
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
	input: 
		sam = "ReadMapping/ema/{sample}/{sample}.nobarcode.sam",
		genome = genomefile
	output: temp("ReadMapping/ema/{sample}/{sample}.nobarcode.bam.tmp")
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
	input: "ReadMapping/ema/{sample}/{sample}.nobarcode.bam.tmp"
	output: 
		bam = temp("ReadMapping/ema/{sample}/{sample}.nobarcode.bam"),
		bai = temp("ReadMapping/ema/{sample}/{sample}.nobarcode.bam.bai")
	log: 
		mdlog = "ReadMapping/ema/log/{sample}.markdup.nobarcode.log",
		stats = "ReadMapping/ema/stats/{sample}.nobarcode.stats",
		flagstat = "ReadMapping/ema/flagstat/{sample}.nobarcode.flagstat"
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
		aln_barcoded = expand("ReadMapping/ema/{{sample}}/{{sample}}.{bin}.bam", bin = ["%03d" % i for i in range(nbins)]),
	output: 
		bam = temp("ReadMapping/ema/{sample}/{sample}.barcoded.bam"),
		bai = temp("ReadMapping/ema/{sample}/{sample}.barcoded.bam.bai")
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message: "Merging barcoded alignments: {wildcards.sample}"
	benchmark: "Benchmark/Mapping/ema/merge.{sample}.txt"
	threads: 10
	shell:
		"""
		sambamba merge -t {threads} {output.bam} {input}
		"""

rule mergestats:
	input: 
		bam = "ReadMapping/ema/{sample}/{sample}.barcoded.bam",
		bai = "ReadMapping/ema/{sample}/{sample}.barcoded.bam.bai"
	log:
		stats = "ReadMapping/ema/stats/{sample}.barcoded.stats",
		flagstat = "ReadMapping/ema/flagstat/{sample}.barcoded.flagstat"
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
	input: "ReadMapping/ema/{sample}/{sample}.barcoded.bam"
	output: 
		filt = temp("ReadMapping/bedfiles/{sample}.bx.bed"),
		unfilt = temp("ReadMapping/bedfiles/{sample}.all.bed")
	message: "Converting to BED format: {wildcards.sample}"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	params: "ReadMapping/ema/{wildcards.sample}/{wildcards.sample}.barcoded.all.bed"
	threads: 1
	shell:
		"""
		utilities/writeBED.pl {input}
		mv {params} {output.unfilt}
		awk '!($4~/A00|B00|C00|D00/)' {output.unfilt} > {output.filt}
		"""

rule BX_stats:
	input: "ReadMapping/bedfiles/{sample}.bx.bed"
	output:	
		molsize = "ReadMapping/bxstats/moleculesize/{sample}.molsize",
		molhist = "ReadMapping/bxstats/moleculesize/{sample}.molsize.hist",
		readsper = "ReadMapping/bxstats/readsperbx/{sample}.readsperbx"
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
		aln_barcoded = "ReadMapping/ema/{sample}/{sample}.barcoded.bam",
		aln_nobarcode = "ReadMapping/ema/{sample}/{sample}.nobarcode.bam",
		idx_barcoded = "ReadMapping/ema/{sample}/{sample}.barcoded.bam.bai",
		idx_nobarcode = "ReadMapping/ema/{sample}/{sample}.nobarcode.bam.bai"
	output: 
		bam = temp("ReadMapping/ema/{sample}.unsort.bam"),
		bai = temp("ReadMapping/ema/{sample}.unsort.bam.bai")
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
		bam = "ReadMapping/ema/{sample}.unsort.bam",
		genome = genomefile
	output: "ReadMapping/ema/{sample}.bam"
	message: "Sorting merged barcoded alignments: {wildcards.sample}"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	threads: 2
	shell:
		"""
		samtools sort -@ {threads} -O bam --reference {input.genome} -l 0 -m 4G -o {output} {input.bam} 2> /dev/null
		"""

rule genome_coords:
	input: genomefile + ".fai"
	output: genomefile + ".bed"
	message: "Creating BED file of genomic coordinates"
	threads: 1
	shell:
		"""
		awk 'BEGIN{{FS=OFS="\\t"}} {{print $1,$2}}' {input} > {output}
		"""

rule bamtobed:
	input: "ReadMapping/ema/{sample}.bam"
	output: temp("ReadMapping/bedfiles/{sample}.total.bed")
	message: "Creating bedfile of all alignments: {wildcards.sample}"
	benchmark: "Benchmark/Mapping/ema/fullbed.{sample}.txt"
	shell:
		"bedtools bamtobed {input}"

rule genome_BX_coverage:
	input:
		geno = genomefile + ".bed",
		bx = "ReadMapping/bedfiles/{sample}.bx.bed",
		alntot = "ReadMapping/bedfiles/{sample}.total.bed"
	output: 
		bx = "ReadMapping/ema/coverage/{sample}.bx.gencov",
		alntot = "ReadMapping/ema/coverage/{sample}.all.gencov"
	message: "Calculating genomic coverage of alignments: {wildcards.sample}"
	shell:
		"""
		bedtools genomecov -i {input.bx}   -g {input.geno} > {output.bx}
		bedtools genomecov -i {input.alntot}   -g {input.geno} > {output.alntot}
		"""

rule index_alignments:
	input: "ReadMapping/ema/{sample}.bam"
	output: "ReadMapping/ema/{sample}.bam.bai"
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
		bam = "ReadMapping/ema/{sample}.bam",
		bai = "ReadMapping/ema/{sample}.bam.bai"
	output:
		stats = report("ReadMapping/ema/stats/{sample}.stats", category="{sample}", subcategory="All Aligments", labels={"Metric": "stats"}),
		flagstat = report("ReadMapping/ema/flagstat/{sample}.flagstat", category="{sample}", subcategory="All Aligments", labels={"Metric": "flagstat"})
	message: "Calculating alignment stats: {wildcards.sample}"
	benchmark: "Benchmark/Mapping/ema/Mergedstats.{sample}.txt"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	shell:
		"""
		samtools stats {input.bam} > {output.stats}
		samtools flagstat {input.bam} > {output.flagstat}
		"""