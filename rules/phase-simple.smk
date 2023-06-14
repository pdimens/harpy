bam_dir           = config["seq_directory"]
samplenames       = config["samplenames"]
variantfile       = config["variantfile"]
pruning           = config["prune"]
molecule_distance = config["molecule_distance"]
extra             = config.get("extra", "") 
fragfile          = "Phase/extractHairs/{sample}.unlinked.frags" if config["noBX"] else "Phase/linkFragments/{sample}.linked.frags"


rule splitbysample:
	input: 
		vcf = variantfile,
		bam = bam_dir + "/{sample}.bam"
	output:
		temp("Phase/input/{sample}.bcf")
	message:
		"Extracting variants: {wildcards.sample}"
	benchmark:
		"Benchmark/Phase/split.{sample}.txt"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	threads: 1
	shell:
		"""
		if grep -q "INFO_SCORE" <(bcftools head {input.vcf}); then
			bcftools view -s {wildcards.sample} -i 'INFO/INFO_SCORE >= 0.2' {input.vcf} 
		else
			bcftools view -s {wildcards.sample} {input.vcf} > {output}
		fi
		#awk '/^#/;/CHROM/ {{OFS="\\t"}}; !/^#/ &&  $10~/^0\\/0/ {{$10="0|0:"substr($10,5);print $0}}; !/^#/ && $10~/^0\\/1/; !/^#/ &&  $10~/^1\\/1/ {{$10="1|1:"substr($10,5);print $0}}; !/^#/ {{print $0}}' > {output}
		"""

rule extractHairs:
	input:
		vcf = "Phase/input/{sample}.bcf",
		bam = bam_dir + "/{sample}.bam"
	output:
		"Phase/extractHairs/{sample}.unlinked.frags"
	log:
		"Phase/extractHairs/logs/{sample}.unlinked.log"
	message:
		"Converting to compact fragment format: {wildcards.sample}"
	wildcard_constraints:
		sample = "[a-zA-Z0-9_-]*"
	benchmark:
		"Benchmark/Phase/extracthairs.{sample}.txt"
	threads: 1
	shell:
		"extractHAIRS --10X 1 --nf 1 --bam {input.bam} --VCF {input.vcf} --out {output} 2> {log}"

rule linkFragments:
	input: 
		bam       = bam_dir + "/{sample}.bam",
		vcf       = "Phase/input/{sample}.bcf",
		fragments = "Phase/extractHairs/{sample}.unlinked.frags"
	output:
		"Phase/linkFragments/{sample}.linked.frags"
	log:
		"Phase/linkFragments/logs/{sample}.linked.log"
	message:
		"Linking fragments: {wildcards.sample}"
	benchmark:
		"Benchmark/Phase/linkfrag.{sample}.txt"
	params:
		d = molecule_distance
	shell:
		"LinkFragments.py  --bam {input.bam} --VCF {input.vcf} --fragments {input.fragments} --out {output} -d {params} > {log} 2>&1"


rule phaseBlocks:
	input:
		vcf       = "Phase/input/{sample}.bcf",
		fragments = fragfile
	output: 
		blocks    = "Phase/phaseBlocks/{sample}.blocks",
		vcf       = "Phase/phaseBlocks/{sample}.blocks.phased.VCF"
	message:
		"Creating phased haplotype blocks: {wildcards.sample}"
	benchmark:
		"Benchmark/Phase/phase.{sample}.txt"
	log:
		"Phase/phaseBlocks/logs/{sample}.blocks.phased.log"
	params: 
		prune = f"--threshold {pruning}" if pruning > 0 else "--no_prune 1",
		extra = extra
	threads: 1
	shell:
		"HAPCUT2 --fragments {input.fragments} --vcf {input.vcf} {params} --out {output.blocks} --nf 1 {params} --error_analysis_mode 1 --call_homozygous 1 --outvcf 1 2> {log}"


rule mergeSamples:
	input: 
		vcf = expand("Phase/phaseBlocks/{sample}.blocks.phased.VCF", sample = samplenames)
	output:
		"Phase/variants.phased.bcf"
	message:
		"Combinging samples into a single BCF file"
	benchmark:
		"Benchmark/Phase/mergesamples.txt"
	threads: 30
	shell:
		"bcftools merge --threads {threads} --output-type b {input.vcf} > {output}"

rule indexFinal:
	input:
		"Phase/variants.phased.bcf"
	output:
		"Phase/variants.phased.bcf.csi"
	benchmark:
		"Benchmark/Phase/finalindex.txt"
	message:
		"Phasing is complete!"
	default_target: True
	shell: 
		"bcftools index {input}"
