bam_dir = config["seq_directory"]
genomefile = config["genomefile"]
samplenames = config["samplenames"] 
extra = config.get("extra", "") 
groupfile = config["groupings"]

bn = os.path.basename(genomefile)
shell("mkdir -p Assembly")
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
		print(f"ERROR: {len(absent)} samples in {infile} not found in {dirn} directory:")
		print(", ".join(absent))
		exit(1)
	return d

popdict = pop_manifest(groupfile, bam_dir, samplenames)
populations = popdict.keys()

rule bamlist:
	output: expand("Variants/leviathan/input/{pop}.list", pop = populations)
	message: "Creating file lists for each population."
	run:
		for p in populations:
			with open(f"Variants/leviathan/input/{p}.list", "w") as fout:
				bamlist = popdict[p]
				for bamfile in bamlist:
					_ = fout.write(bamfile + "\n")

rule merge_alignments:
	input: 
		bamlist = "Variants/leviathan/input/{population}.list",
		bamfiles = lambda wc: expand("{sample}", sample = popdict[wc.population]) 
	output: temp("Variants/leviathan/input/{population}.bam")
	message: "Merging alignments: Population {wildcards.population}"
	shell:
		"""
		samtools merge -b {input} -o {output}        
		"""

rule index_merged:
	input: "Variants/leviathan/input/{population}.bam"
	output: temp("Variants/leviathan/input/{population}.bam.bai")
	message: "Indexing alignments: Population {wildcards.population}"
	wildcard_constraints:
		population = "[a-zA-Z0-9_-]*"
	shell:
		"""
		sambamba index {input} {output}
		"""

rule keep_validBX:
	input: "Variants/leviathan/input/{population}.bam"
	output: "Variants/leviathan/input/{population}.bx.valid.bam"
	message: "Keeping only alignments with valid BX barcodes: {wildcards.population}"
	wildcard_constraints:
		population = "[a-zA-Z0-9_-]*"
	shell:
		"""
		utilities/filterBXBAM.py --valid --input {input}
		"""

rule index_valid:
	input: "Variants/leviathan/input/{population}.bx.valid.bam"
	output: "Variants/leviathan/input/{population}.bx.valid.bam.bai"
	message: "Indexing barcodes: Population {wildcards.population}"
	benchmark: "Benchmark/Variants/leviathan/indexbam.{population}.txt"
	wildcard_constraints:
		population = "[a-zA-Z0-9_-]*"
	threads: 1
	shell:
		"""
		sambamba index {input} {output}
		"""

rule index_barcode:
	input: 
		bam = "Variants/leviathan/input/{population}.bx.valid.bam",
		bai = "Variants/leviathan/input/{population}.bx.valid.bam.bai"
	output: temp("Variants/leviathan/lrezIndexed/{population}.bci")
	message: "Indexing barcodes: Population {wildcards.population}"
	benchmark: "Benchmark/Variants/leviathan/indexbc.{population}.txt"
	threads: 4
	shell:
		"""
		LRez index bam -p -b {input.bam} -o {output} --threads {threads}
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

rule leviathan_variantcall:
	input:
		bam = "Variants/leviathan/input/{population}.bx.valid.bam",
		bai = "Variants/leviathan/input/{population}.bx.valid.bam.bai",
		bc_idx = "Variants/leviathan/lrezIndexed/{population}.bci",
		genome = f"Assembly/{genomefile}"
	output: vcf = pipe("Variants/leviathan/{population}.vcf")
	log:  
		runlog = "Variants/leviathan/logs/{population}.leviathan.log",
		candidates = "Variants/leviathan/logs/{population}.candidates"
	message: "Calling variants: Population {wildcards.population}"
	benchmark: "Benchmark/Variants/leviathan/variantcall.{population}.txt"
	params:
		extra = extra
	threads: 3
	shell:
		"""
		LEVIATHAN -b {input.bam} -i {input.bc_idx} {params} -g {input.genome} -o {output} -t {threads} --candidates {log.candidates} 2> {log.runlog}
		"""

rule sort_bcf:
	input: "Variants/leviathan/{population}.vcf"
	output: "Variants/leviathan/{population}.bcf"
	message: "Sorting and converting to BCF: Population {wildcards.population}"
	threads: 1
	params: "{wildcards.population}"
	benchmark: "Benchmark/Variants/leviathan/sortbcf.{population}.txt"
	shell:        
		"""
		bcftools sort -Ob --output {output} {input} 2> /dev/null
		"""

rule index_bcf:
	input: "Variants/leviathan/{population}.bcf"
	output: "Variants/leviathan/{population}.bcf.csi"
	message: "Indexing: Population {input}"
	benchmark: "Benchmark/Variants/leviathan/indexbcf.{population}.txt"
	threads: 1
	shell:
		"""
		bcftools index --output {output} {input}
		"""

rule sv_stats:
	input: 
		bcf = "Variants/leviathan/{population}.bcf",
		idx = "Variants/leviathan/{population}.bcf.csi"
	output: "Variants/leviathan/stats/{population}.sv.stats"
	message: "Getting stats: Population {input.bcf}"
	benchmark: "Benchmark/Variants/leviathan/stats.{population}.txt"
	threads: 1
	shell:
		"""
		bcftools stats {input.bcf} > {output}
		"""

rule all_bcfs:
	input: 
		bcf = expand("Variants/leviathan/{pop}.bcf", pop = populations),
		stats = expand("Variants/leviathan/stats/{pop}.sv.stats", pop = populations)
	message: "Variant calling is complete!"
	default_target: True