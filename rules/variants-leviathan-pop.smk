import sys
import os

bam_dir = config["seq_directory"]
genomefile = config["genomefile"]
samplenames = config["samplenames"] 
extra = config.get("extra", "") 
groupfile = config["groupings"]

bn = os.path.basename(genomefile)
os.makedirs("Assembly", exist_ok = True)
if not os.path.exists(f"Assembly/{bn}"):
	shell(f"ln -sr {genomefile} Assembly/{bn}")

# create dictionary of population => filenames
## this makes it easier to set the snakemake rules/wildcards
## exits with an error if the groupfile has samples not in the bam folder
def pop_manifest(infile, dirn, sampnames):
	d = dict()
	absent = []
	with open(infile) as f:
		for line in f:
			samp, pop = line.rstrip().split()
			if samp not in sampnames:
				absent.append(samp)
			samp = f"{dirn}/{samp}.bam"
			if pop not in d.keys():
				d[pop] = [samp]
			else:
				d[pop].append(samp)
	if absent:
		sys.tracebacklimit = 0
		raise ValueError(f"{len(absent)} sample(s) in \033[1m{infile}\033[0m not found in \033[1m{dirn}\033[0m directory:\n\033[33m" + ", ".join(absent) + "\033[0m")
	return d

popdict = pop_manifest(groupfile, bam_dir, samplenames)
populations = popdict.keys()

rule bamlist:
	output:
		expand("Variants/leviathan-pop/input/{pop}.list", pop = populations)
	message:
		"Creating file lists for each population."
	run:
		for p in populations:
			with open(f"Variants/leviathan-pop/input/{p}.list", "w") as fout:
				bamlist = popdict[p]
				for bamfile in bamlist:
					_ = fout.write(bamfile + "\n")

rule merge_populations:
	input: 
		bamlist = "Variants/leviathan-pop/input/{population}.list",
		bamfiles = lambda wc: expand("{sample}", sample = popdict[wc.population]) 
	output:
		temp("Variants/leviathan-pop/input/{population}.bam")
	message:
		"Merging alignments: Population {wildcards.population}"
	shell:
		"samtools merge -b {input} -o {output}"

rule index_merged:
	input:
		"Variants/leviathan-pop/input/{population}.bam"
	output:
		temp("Variants/leviathan-pop/input/{population}.bam.bai")
	message:
		"Indexing merged alignments: Population {wildcards.population}"
	wildcard_constraints:
		population = "[a-zA-Z0-9_-]*"
	shell:
		"sambamba index {input} {output} 2> /dev/null"

rule index_barcode:
	input: 
		bam = "Variants/leviathan-pop/input/{population}.bam",
		bai = "Variants/leviathan-pop/input/{population}.bam.bai"
	output:
		temp("Variants/leviathan-pop/lrezIndexed/{population}.bci")
	message:
		"Indexing barcodes: Population {wildcards.population}"
	benchmark:
		"Benchmark/Variants/leviathan-pop/indexbc.{population}.txt"
	threads: 4
	shell:
		"LRez index bam -p -b {input.bam} -o {output} --threads {threads}"

rule index_genome:
	input:
		genomefile
	output: 
		asm = f"Assembly/{genomefile}",
		idx = multiext(f"Assembly/{genomefile}", ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
	message:
		"Indexing {input}"
	log:
		f"Assembly/{genomefile}.idx.log"
	shell: 
		"""
		ln -sr {input} {output.asm}
		bwa index {output.asm} 2> {log}
		samtools faidx --fai-idx {output.asm}.fai {output.asm} 2>> {log}
		"""

rule leviathan_variantcall:
	input:
		bam = "Variants/leviathan-pop/input/{population}.bam",
		bai = "Variants/leviathan-pop/input/{population}.bam.bai",
		bc_idx = "Variants/leviathan-pop/lrezIndexed/{population}.bci",
		genome = f"Assembly/{genomefile}"
	output:
		pipe("Variants/leviathan-pop/{population}.vcf")
	log:  
		runlog = "Variants/leviathan-pop/logs/{population}.leviathan.log",
		candidates = "Variants/leviathan-pop/logs/{population}.candidates"
	message:
		"Calling variants: Population {wildcards.population}"
	benchmark:
		"Benchmark/Variants/leviathan-pop/variantcall.{population}.txt"
	params:
		extra = extra
	threads: 3
	shell:
		"LEVIATHAN -b {input.bam} -i {input.bc_idx} {params} -g {input.genome} -o {output} -t {threads} --candidates {log.candidates} 2> {log.runlog}"

rule sort_bcf:
	input:
		"Variants/leviathan-pop/{population}.vcf"
	output:
		"Variants/leviathan-pop/{population}.bcf"
	message:
		"Sorting and converting to BCF: Population {wildcards.population}"
	threads: 1
	params:
		"{wildcards.population}"
	benchmark:
		"Benchmark/Variants/leviathan-pop/sortbcf.{population}.txt"
	shell:        
		"bcftools sort -Ob --output {output} {input} 2> /dev/null"

rule sv_stats:
	input: 
		"Variants/leviathan-pop/{population}.bcf"
	output:
		"Variants/leviathan-pop/reports/stats/{population}.sv.stats"
	message:
		"Getting stats: Population {input}"
	benchmark:
		"Benchmark/Variants/leviathan-pop/stats.{population}.txt"
	threads: 1
	shell:
		"""
		echo -e "population\\tcontig\\tposition_start\\tposition_end\\tlength\\ttype\\tn_barcodes\\tn_pairs" > {output}
		bcftools query -f '{wildcards.population}\\t%CHROM\\t%POS\\t%END\\t%SVLEN\\t%SVTYPE\\t%BARCODES\\t%PAIRS\\n' {input} >> {output}
		"""

rule sv_report_bypop:
	input:	
		statsfile = "Variants/leviathan-pop/reports/stats/{population}.sv.stats",
		faidx = f"Assembly/{genomefile}.fai"
	output:
		"Variants/leviathan-pop/reports/{population}.sv.html"
	message:
		"Generating SV report for all populations"
	script:
		"reportLeviathan.Rmd"


rule sv_report:
	input:	
		statsfiles = expand("Variants/leviathan-pop/reports/stats/{pop}.sv.stats", pop = populations),
		faidx = f"Assembly/{genomefile}.fai"
	output:
		"Variants/leviathan-pop/reports/SV.summary.html"
	message:
		"Generating SV report for all populations"
	script:
		"reportLeviathanPop.Rmd"

rule all_bcfs:
	input: 
		bcf = expand("Variants/leviathan-pop/{pop}.bcf", pop = populations),
		stats = expand("Variants/leviathan-pop/reports/stats/{pop}.sv.stats", pop = populations),
		popreports = expand("Variants/leviathan-pop/reports/{pop}.sv.html", pop = populations),
		finalreport = "Variants/leviathan-pop/reports/SV.summary.html"
	default_target: True
	message:
		"Variant calling is complete!"