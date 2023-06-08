import sys
import os

bam_dir 	= config["seq_directory"]
genomefile 	= config["genomefile"]
samplenames = config["samplenames"] 
extra 		= config.get("extra", "") 
groupfile 	= config["groupings"]
bn 			= os.path.basename(genomefile)
outdir      = "Variants/leviathan-pop"

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
		raise ValueError(f"{len(absent)} sample(s) in \033[1m{infile}\033[0m not found in \033[1m{dirn}\033[0m directory:\n\033[33m" + ", ".join(absent) + "\033[0m" + "\n")
	return d

popdict = pop_manifest(groupfile, bam_dir, samplenames)
populations = popdict.keys()

rule bamlist:
	output:
		expand(outdir + "/input/{pop}.list", pop = populations)
	message:
		"Creating file lists for each population."
	run:
		for p in populations:
			bamlist = popdict[p]
			with open(f"{outdir}/input/{p}.list", "w") as fout:
				for bamfile in bamlist:
					_ = fout.write(bamfile + "\n")

rule merge_populations:
	input: 
		bamlist  = outdir + "/input/{population}.list",
		bamfiles = lambda wc: expand("{sample}", sample = popdict[wc.population]) 
	output:
		temp(outdir + "/input/{population}.bam")
	message:
		"Merging alignments: Population {wildcards.population}"
	shell:
		"samtools merge -b {input} -o {output}"

rule index_merged:
	input:
		outdir + "/input/{population}.bam"
	output:
		temp(outdir + "/input/{population}.bam.bai")
	message:
		"Indexing merged alignments: Population {wildcards.population}"
	wildcard_constraints:
		population = "[a-zA-Z0-9_-]*"
	shell:
		"sambamba index {input} {output} 2> /dev/null"

rule index_barcode:
	input: 
		bam = outdir + "/input/{population}.bam",
		bai = outdir + "/input/{population}.bam.bai"
	output:
		temp(outdir + "/lrezIndexed/{population}.bci")
	message:
		"Indexing barcodes: Population {wildcards.population}"
	benchmark:
		"Benchmark/Variants/leviathan-pop/indexbc.{population}.txt"
	threads: 4
	shell:
		"LRez index bam -p -b {input.bam} -o {output} --threads {threads}"

rule link_genome:
	input:
		genomefile
	output: 
		f"Assembly/{bn}"
	message:
		"Symlinking {input} to Assembly/"
	shell: 
		"ln -sr {input} {output}"

rule index_faidx_genome:
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

rule leviathan_variantcall:
	input:
		bam    = outdir + "/input/{population}.bam",
		bai    = outdir + "/input/{population}.bam.bai",
		bc_idx = outdir + "/lrezIndexed/{population}.bci",
		genome = f"Assembly/{bn}"
	output:
		pipe(outdir + "/{population}.vcf")
	log:  
		runlog     = outdir + "/logs/{population}.leviathan.log",
		candidates = outdir + "/logs/{population}.candidates"
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
		outdir + "/{population}.vcf"
	output:
		outdir + "/{population}.bcf"
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
		outdir + "/{population}.bcf"
	output:
		outdir + "/reports/stats/{population}.sv.stats"
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
		statsfile = outdir + "/reports/stats/{population}.sv.stats",
		bcf       = outdir + "/{population}.bcf"
	output:
		outdir + "/reports/{population}.sv.html"
	message:
		"Generating SV report: population {wildcards.population}"
	script:
		"reportLeviathan.Rmd"


rule sv_report:
	input:	
		faidx      = f"Assembly/{bn}.fai",
		statsfiles = expand(outdir + "/reports/stats/{pop}.sv.stats", pop = populations)
	output:
		outdir + "/reports/leviathan.pop.summary.html"
	message:
		"Generating SV report for all populations"
	script:
		"reportLeviathanPop.Rmd"

rule all_bcfs:
	input: 
		bcf       = expand(outdir + "/{pop}.bcf", pop = populations),
		stats     = expand(outdir + "/reports/stats/{pop}.sv.stats", pop = populations),
		popreport = expand(outdir + "/reports/{pop}.sv.html", pop = populations),
		report    = outdir + "/reports/leviathan.pop.summary.html"
	default_target: True
	message:
		"Variant calling is complete!"