import os

seq_dir 	= config["seq_directory"]
nbins 		= config["EMA_bins"]
genomefile 	= config["genomefile"]
Rsep 		= config["Rsep"]
fqext 		= config["fqext"]
samplenames = config["samplenames"]
extra 		= config.get("extra", "") 
bn 			= os.path.basename(genomefile)
outdir      = "Align/ema"

rule create_reports:
	input: 
		expand(outdir + "/{sample}.bam", sample = samplenames),
		expand(outdir + "/{sample}.bam.bai", sample = samplenames),
		expand(outdir + "/stats/moleculesize/{sample}.{ext}", sample = samplenames, ext = ["molsize.gz", "molsize.hist"]),
		expand(outdir + "/stats/readsperbx/{sample}.readsperbx", sample = samplenames),
		expand(outdir + "/stats/coverage/{sample}.cov.html", sample = samplenames),
		outdir + "/stats/reads.bxstats.html",
		outdir + "/stats/samtools_stats/alignment.stats.html",
		outdir + "/stats/samtools_flagstat/alignment.flagstat.html"
	message:
		"Read mapping completed!"
	benchmark:
		"Benchmark/Mapping/ema/report.txt"
	default_target: True

rule link_genome:
	input:
		genomefile
	output: 
		f"Assembly/{bn}"
	message:
		"Symlinking {input} to Assembly/"
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
		counts = outdir + "/count/{sample}.ema-ncnt",
		logs   = temp(outdir + "/count/logs/{sample}.count.log")
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message:
		"Counting barcode frequency: {wildcards.sample}"
	benchmark:
		"Benchmark/Mapping/ema/Count.{sample}.txt"
	params:
		prefix = lambda wc: outdir + "/count/" + wc.get("sample")
	threads: 1
	shell:
		"seqfu interleave -1 {input.forward_reads} -2 {input.reverse_reads} | ema-h count -p -o {params} 2> {output.logs}"

rule beadtag_summary:
	input: 
		countlog = expand(outdir + "/count/logs/{sample}.count.log", sample = samplenames)
	output:
		outdir + "/stats/reads.bxstats.html"
	message:
		"Creating sample barcode validation report"
	benchmark:
		"Benchmark/Mapping/ema/beadtagsummary.txt"
	script:
		"reportBxCount.Rmd"

rule preprocess_ema:
	input: 
		forward_reads = seq_dir + "/{sample}" + f".{Rsep[0]}.{fqext}",
		reverse_reads = seq_dir + "/{sample}" + f".{Rsep[1]}.{fqext}",
		emacounts     = outdir + "/count/{sample}.ema-ncnt"
	output: 
		bins       	  = temp(expand(outdir + "/preproc/{{sample}}/ema-bin-{bin}", bin = ["%03d" % i for i in range(nbins)])),
		unbarcoded    = temp(outdir + "/preproc/{sample}/ema-nobc")
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	log:
		outdir + "/preproc/logs/{sample}.preproc.log"
	message:
		"Preprocessing for EMA mapping: {wildcards.sample}"
	benchmark:
		"Benchmark/Mapping/ema/Preproc.{sample}.txt"
	threads:
		2
	params:
		outdir = lambda wc: outdir + "/preproc/" + wc.get("sample"),
		bins   = nbins
	shell:
		"seqfu interleave -1 {input.forward_reads} -2 {input.reverse_reads} | ema-h preproc -p -n {params.bins} -t {threads} -o {params.outdir} {input.emacounts} 2>&1 | cat - > {log}"

rule align_ema:
	input:
		readbin    = outdir + "/preproc/{sample}/ema-bin-{bin}",
		genome 	   = f"Assembly/{bn}",
		genome_idx = multiext(f"Assembly/{bn}", ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
	output:
		alignment  = temp(outdir + "/align/{sample}/{sample}.{bin}.bam")
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message:
		"Aligning barcoded sequences: {wildcards.sample}-{wildcards.bin}"
	benchmark:
		"Benchmark/Mapping/ema/Align.{sample}.{bin}.txt"
	params: 
		quality = config["quality"],
		extra = extra
	threads: 8
	shell:
		"""
		EMATHREADS=$(( {threads} - 2 ))
		ema-h align -t $EMATHREADS {params.extra} -d -p haptag -r {input.genome} -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" -s {input.readbin} 2> /dev/null |
		samtools view -h -F 4 -q {params.quality} - | 
		samtools sort --reference {input.genome} -O bam -m 4G -o {output} - 2> /dev/null
		"""

rule align_nobarcode:
	input:
		reads      = outdir + "/preproc/{sample}/ema-nobc",
		genome 	   = f"Assembly/{bn}",
		genome_idx = multiext(f"Assembly/{bn}", ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
	output: 
		samfile    = temp(outdir + "/align/{sample}/{sample}.nobarcode.bam.tmp")
	benchmark:
		"Benchmark/Mapping/ema/bwaAlign.{sample}.txt"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	params:
		quality = config["quality"]
	message:
		"Aligning unbarcoded sequences: {wildcards.sample}"
	threads: 8
	shell:
		"""
		BWATHREADS=$(( {threads} - 2 ))
		bwa mem -t $BWATHREADS -C -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" {input.genome} {input.reads} 2> /dev/null |
		samtools view -h -F 4 -q {params.quality} | 
		samtools sort -O bam -m 4G --reference {input.genome} -o {output} 2> /dev/null
		"""

rule markduplicates:
	input:
		bam      = outdir + "/align/{sample}/{sample}.nobarcode.bam.tmp"
	output: 
		bam      = temp(outdir + "/align/{sample}/{sample}.nobarcode.bam"),
		bai      = temp(outdir + "/align/{sample}/{sample}.nobarcode.bam.bai")
	log: 
		mdlog    = outdir + "/stats/markduplicates/{sample}.markdup.nobarcode.log",
		stats    = outdir + "/stats/samtools_stats/{sample}.nobarcode.stats",
		flagstat = outdir + "/stats/samtools_flagstat/{sample}.nobarcode.flagstat"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message:
		"Marking duplicates in unbarcoded alignments: {wildcards.sample}"
	benchmark:
		"Benchmark/Mapping/ema/markdup.{sample}.txt"
	threads: 2
	shell:
		"""
		sambamba markdup -t {threads} -l 4 {input} {output.bam} 2> {log.mdlog}
		samtools stats {output.bam} > {log.stats}
		samtools flagstat {output.bam} > {log.flagstat}
		"""   

rule merge_barcoded:
	input:
		aln_barcoded = expand(outdir + "/align/{{sample}}/{{sample}}.{bin}.bam", bin = ["%03d" % i for i in range(nbins)]),
	output: 
		bam 		 = temp(outdir + "/align/barcoded/{sample}.barcoded.bam"),
		bai 		 = temp(outdir + "/align/barcoded/{sample}.barcoded.bam.bai")
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message:
		"Merging barcoded alignments: {wildcards.sample}"
	benchmark:
		"Benchmark/Mapping/ema/merge.{sample}.txt"
	threads: 10
	shell:
		"sambamba merge -t {threads} -l 4 {output.bam} {input} 2> /dev/null"

#rule secondary2split:
#	input:
#		bam = outdir + "/align/barcoded/{sample}.barcoded.sec.bam",
#		bai = outdir + "/align/barcoded/{sample}.barcoded.sec.bam.bai"
#	output:
#		bam = temp(outdir + "/align/barcoded/{sample}.barcoded.bam"),
#		bai = temp(outdir + "/align/barcoded/{sample}.barcoded.bam.bai")
#	wildcard_constraints:
#		sample = "[a-zA-Z0-9_-]*"
#	message:
#		"Converting Secondary SAM flags to Split flags: {wildcards.sample}"
#	shell:
#		"secondary2split.py {input.bam} {output.bam}"

rule bcstats:
	input: 
		bam      = outdir + "/align/barcoded/{sample}.barcoded.bam",
		bai      = outdir + "/align/barcoded/{sample}.barcoded.bam.bai"
	log:
		stats    = outdir + "/stats/samtools_stats/{sample}.barcoded.stats",
		flagstat = outdir + "/stats/samtools_flagstat/{sample}.barcoded.flagstat"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message:
		"Indexing merged barcoded alignemnts: {wildcards.sample}"
	benchmark:
		"Benchmark/Mapping/ema/indexmerge.{sample}.txt"
	shell:
		"""
		samtools stats {input.bam} > {log.stats}
		samtools flagstat {input.bam} > {log.flagstat}
		"""

rule alignment_coverage:
	input: 
		bed     = f"Assembly/{bn}.bed",
		nobx    = outdir + "/align/{sample}/{sample}.nobarcode.bam",
		nobxbai = outdir + "/align/{sample}/{sample}.nobarcode.bam.bai",
		bx      = outdir + "/align/barcoded/{sample}.barcoded.bam",
		bxbai   = outdir + "/align/barcoded/{sample}.barcoded.bam.bai"
	output: 
		outdir + "/stats/coverage/data/{sample}.cov.gz"
	message:
		"Calculating genomic coverage: {wildcards.sample}"
	threads: 2
	shell:
		"samtools bedcov -c {input.bed} {input.bx} {input.nobx} | gzip > {output}"

rule gencovBX_report:
	input: 
		outdir + "/stats/coverage/data/{sample}.cov.gz",
	output:
		outdir + "/stats/coverage/{sample}.cov.html"
	message:
		"Creating report of alignment coverage: {wildcards.sample}"
	script:
		"reportEmaGencov.Rmd"

rule BEDconvert:
	input:
		outdir + "/align/barcoded/{sample}.barcoded.bam"
	output: 
		unfilt = temp(outdir + "/bedfiles/{sample}.all.bed"),
		bx     = temp(outdir + "/bedfiles/{sample}.bx.bed")
	message:
		"Converting to BED format: {wildcards.sample}"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	params: lambda wc: outdir + "/align/" + wc.get("sample") + "/" + wc.get("sample") + ".all.bed"
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

rule merge_alignments:
	input:
		aln_barcoded  = outdir + "/align/barcoded/{sample}.barcoded.bam",
		idx_barcoded  = outdir + "/align/barcoded/{sample}.barcoded.bam.bai",
		aln_nobarcode = outdir + "/align/{sample}/{sample}.nobarcode.bam",
		idx_nobarcode = outdir + "/align/{sample}/{sample}.nobarcode.bam.bai"
	output: 
		bam 		  = temp(outdir + "/align/{sample}.unsort.bam"),
		bai 		  = temp(outdir + "/align/{sample}.unsort.bam.bai")
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	message:
		"Merging all alignments: {wildcards.sample}"
	benchmark:
		"Benchmark/Mapping/ema/mergebc_nobc.{sample}.txt"
	threads: 10
	shell:
		"sambamba merge -t {threads} {output.bam} {input.aln_barcoded} {input.aln_nobarcode} 2> /dev/null"

rule sort_merge:
	input:
		bam    = outdir + "/align/{sample}.unsort.bam",
		genome = f"Assembly/{bn}"
	output:
		outdir + "/{sample}.bam"
	message:
		"Sorting merged barcoded alignments: {wildcards.sample}"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	threads: 2
	priority: 1
	shell:
		"samtools sort -@ {threads} -O bam --reference {input.genome} -m 4G -o {output} {input.bam} 2> /dev/null"

rule index_alignments:
	input: 
		outdir + "/{sample}.bam"
	output:
		outdir + "/{sample}.bam.bai"
	message:
		"Indexing: {input}"
	benchmark:
		"Benchmark/Mapping/ema/IndexMerged.{sample}.txt"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	shell:
		"sambamba index {input} {output} 2> /dev/null"

rule alignment_stats:
	input: 		
		bam      = outdir + "/{sample}.bam",
		bai      = outdir + "/{sample}.bam.bai"
	output:
		stats    = outdir + "/stats/samtools_stats/{sample}.stats",
		flagstat = outdir + "/stats/samtools_flagstat/{sample}.flagstat"
	message:
		"Calculating alignment stats: {wildcards.sample}"
	benchmark:
		"Benchmark/Mapping/ema/Mergedstats.{sample}.txt"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	shell:
		"""
		samtools stats {input.bam} > {output.stats}
		samtools flagstat {input.bam} > {output.flagstat}
		"""

rule samtools_reports:
	input: 
		expand(outdir + "/stats/samtools_{ext}/{sample}.{ext}", sample = samplenames, ext = ["stats", "flagstat"]),
	output: 
		stats    = outdir + "/stats/samtools_stats/alignment.stats.html",
		flagstat = outdir + "/stats/samtools_flagstat/alignment.flagstat.html"
	message:
		"Summarizing samtools stats and flagstats"
	benchmark:
		"Benchmark/Mapping/ema/report.txt"
	shell:
		"""
		multiqc Align/ema/stats/samtools_stats    --force --quiet --no-data-dir --filename {output.stats} 2> /dev/null
		multiqc Align/ema/stats/samtools_flagstat --force --quiet --no-data-dir --filename {output.flagstat} 2> /dev/null
		"""