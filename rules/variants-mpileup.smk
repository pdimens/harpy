import os
import sys

bam_dir 	= config["seq_directory"]
genomefile 	= config["genomefile"]
bn          = os.path.basename(genomefile)
groupings 	= config.get("groupings", None)
ploidy 		= config["ploidy"]
samplenames = config["samplenames"]
mp_extra 	= config.get("extra", "") 
outdir      = "Variants/mpileup"

if groupings is not None:
	absent = []
	with open(groupings) as f:
		for line in f:
			samp, pop = line.rstrip().split()
			if samp not in samplenames:
				absent.append(samp)
	if absent:
		sys.tracebacklimit = 0
		raise ValueError(f"{len(absent)} sample(s) in \033[1m{groupings}\033[0m not found in \033[1m{bam_dir}\033[0m directory:\n\033[33m" + ", ".join(absent) + "\033[0m")

def faidx_contignames(infile):
	bn = os.path.basename(infile)
	os.makedirs("Assembly", exist_ok = True)
	if not os.path.exists(f"Assembly/{bn}"):
		shell(f"ln -sr {infile} Assembly/{bn}")
	if not os.path.exists(f"Assembly/{bn}.fai"):
		print(f"Assembly/{bn}.fai not found, indexing {bn} with samtools faidx", file = sys.stderr)
		subprocess.run(["samtools","faidx", "--fai-idx", f"Assembly/{bn}.fai", infile, "2>", "/dev/null"])
	with open(f"Assembly/{bn}.fai") as f:
		lines = [line.rstrip().split("\t")[0] for line in f]
	return lines

contigs   = faidx_contignames(genomefile)
dict_cont = dict(zip(contigs, contigs))

rule index_alignments:
	input:
		bam_dir + "/{sample}.bam"
	output:
		bam_dir + "/{sample}.bam.bai"
	message:
		"Indexing alignments: {wildcards.sample}"
	benchmark:
		"Benchmark/Variants/mpileup/indexbam.{sample}.txt"
	shell:
		"sambamba index {input} {output} 2> /dev/null"

rule bam_list:
	input: 
		bam = expand(bam_dir + "/{sample}.bam", sample = samplenames),
		bai = expand(bam_dir + "/{sample}.bam.bai", sample = samplenames)
	output:
		outdir + "/logs/samples.files"
	message:
		"Creating list of alignment files"
	benchmark:
		"Benchmark/Variants/mpileup/bamlist.txt"
	run:
		with open(output[0], "w") as fout:
			for bamfile in input.bam:
				_ = fout.write(bamfile + "\n")

rule samplenames:
	output:
		outdir + "/logs/samples.names"
	message:
		"Creating list of sample names"
	run:
		with open(output[0], "w") as fout:
			for samplename in samplenames:
				_ = fout.write(samplename + "\n")		

rule mpileup:
	input:
		bamlist = outdir + "/logs/samples.files",
        genome  = f"Assembly/{bn}"
	output: 
		pipe(outdir + "/{part}.mp.bcf")
	params: 
		lambda wc: dict_cont[wc.part]
	message: 
		"Finding variants: {wildcards.part}"
	log: 
		outdir + "/logs/{part}.mpileup.log"
	benchmark: 
		"Benchmark/Variants/mpileup/mpileup.{part}.txt"
	params:
		region = "{wildcards.part}",
		extra = mp_extra
	shell:
		"bcftools mpileup --fasta-ref {input.genome} --region {params} --bam-list {input.bamlist} --annotate AD --output-type b > {output} 2> {log}"

rule call_genotypes:
	input:
		outdir + "/{part}.mp.bcf"
	output:
		temp(outdir + "/call/{part}.bcf")
	message:
		"Calling genotypes: {wildcards.part}"
	benchmark:
		"Benchmark/Variants/mpileup/call.{part}.txt"
	log:
		outdir + "/logs/{part}.call.log"
	threads: 2
	params: 
		groupsamples = '' if groupings is None else f"--group-samples {groupings}",
		ploidy = f"--ploidy {ploidy}"
	shell:
		"bcftools call --multiallelic-caller {params} --variants-only --output-type b {input} | bcftools sort - --output {output} 2> /dev/null"

rule index_bcf:
	input: 
		bcf     = outdir + "/call/{part}.bcf",
		samples = outdir + "/logs/samples.names",
		genome  = f"Assembly/{genomefile}"
	output:
		temp(outdir + "/call/{part}.bcf.csi")
	log:
		outdir + "/stats/{part}.stats"
	message:
		"Indexing: {wildcards.part}"
	benchmark:
		"Benchmark/Variants/mpileup/indexbcf.{part}.txt"
	threads: 4
	shell:
		"""
		bcftools index --threads {threads} --output {output} {input.bcf}
		bcftools stats -S {input.samples} --fasta-ref {input.genome} {input.bcf} > {log}
		"""

rule combine_bcfs:
	input: 
		bcf     = expand(outdir + "/call/{part}.bcf", part = contigs),
		idx     = expand(outdir + "/call/{part}.bcf.csi", part = contigs),
		genome  = f"Assembly/{genomefile}",
		samples = outdir + "/logs/samples.names"
	output: 
		bcf     = outdir + "/variants.raw.bcf",
		idx     = outdir + "/variants.raw.bcf.csi",
		stats   = outdir + "/stats/variants.raw.stats"
	message:
		"Merging all BCFs into: {output.bcf}"
	benchmark:
		"Benchmark/Variants/mpileup/merge.txt"
	threads: 50
	shell:
		"""
		bcftools concat --threads {threads} --output-type b --naive {input.bcf} > {output.bcf} 2> /dev/null
		bcftools index --output {output.idx} {output.bcf}
		bcftools stats -S {input.samples} --fasta-ref {input.genome} {output.bcf} > {output.stats}
		"""

rule normalize_bcf:
	input: 
		genome  = f"Assembly/{genomefile}",
		bcf     = outdir + "/variants.raw.bcf",
		samples = outdir + "/logs/samples.names"
	output:
		bcf     = outdir + "/variants.normalized.bcf",
		idx     = outdir + "/variants.normalized.bcf.csi",
		stats   = outdir + "/stats/variants.normalized.stats"
	message: 
		"Normalizing the called variants"
	threads: 2
	shell:
		"""
		bcftools norm -d none -f {input.genome} {input.bcf} | bcftools norm -m -any -N -Ob > {output.bcf}
		bcftools index --output {output.idx} {output.bcf}
		bcftools stats -S {input.samples} --fasta-ref {input.genome} {output.bcf} > {output.stats}
		"""

rule bcfreport:
	input:
		outdir + "/stats/variants.raw.stats"
	output:
		outdir + "/stats/variants.raw.html"
	message:
		"Generating bcftools report: variants.raw.bcf"
	benchmark:
		"Benchmark/Variants/mpileup/reports.txt"
	script:
		"reportBcftools.Rmd"

rule bcfreportnorm:
	input:
		outdir + "/stats/variants.normalized.stats"
	output:
		outdir + "/stats/variants.normalized.html"
	message:
		"Generating bcftools report: variants.normalized.bcf"
	script:
		"reportBcftools.Rmd"

rule all:
	input: 
		expand(outdir + "/variants.{file}.bcf", file = ["raw","normalized"]),
		expand(outdir + "/stats/variants.{file}.html", file = ["raw","normalized"])
	message:
		"Variant calling is complete!"
	default_target: True

#with open(f"{outdir}/logs/variants.params", "w") as f:
#	_ = f.write("The harpy variants module ran using these parameters:\n\n")
#	_ = f.write("bcftools mpileup --fasta-ref GENOME --region CONTIG " + mp_extra + " --bam-list BAMS --annotate AD --output-type b\n")
#	gp = '' if groupings is None else f"--group-samples {groupings} " + f"--ploidy {ploidy}"
#	_ = f.write("bcftools call --multiallelic-caller " + gp + " --variants-only --output-type b | bcftools sort -\n")
#	_ = f.write("bcftools concat --output-type b --naive\n")
#	_ = f.write("bcftools norm -d none | bcftools norm -m -any -N -Ob")