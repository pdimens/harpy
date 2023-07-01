import os
import sys

bam_dir 	= config["seq_directory"]
genomefile 	= config["genomefile"]
bn          = os.path.basename(genomefile)
groupings 	= config.get("groupings", None)
ploidy 		= config["ploidy"]
samplenames = config["samplenames"]
extra 	    = config.get("extra", "") 
outdir      = "Variants/freebayes"
chunksize   = config["windowsize"]

# create a python list of regions instead of creating a multitude of files
def createregions(infile, window):
    bn = os.path.basename(infile)
    os.makedirs("Assembly", exist_ok = True)
    if not os.path.exists(f"Assembly/{bn}"):
        shell(f"ln -sr {infile} Assembly/{bn}")
    if not os.path.exists(f"Assembly/{bn}.fai"):
        print(f"Assembly/{bn}.fai not found, indexing {bn} with samtools faidx", file = sys.stderr)
        subprocess.run(["samtools","faidx", "--fai-idx", f"Assembly/{bn}.fai", infile, "2>", "/dev/null"])
    with open(f"Assembly/{bn}.fai") as fai:
        bedregion = []
        while True:
            # Get next line from file
            line = fai.readline()
            # if line is empty, end of file is reached
            if not line:
                break
            # split the line by tabs
            lsplit = line.split()
            contig = lsplit[0]
            c_len = int(lsplit[1])
            start = 0
            end = window
            starts = [0]
            ends = [window]
            while end < c_len:
                end = end + window if (end + window) < c_len else c_len
                ends.append(end)
                start += window
                starts.append(start)
            for (startpos, endpos) in zip (starts,ends):
                bedregion.append(f"{contig}:{startpos}-{endpos}")
        return bedregion

_regions   = createregions(genomefile, chunksize)
regions = dict(zip(_regions, _regions))

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

rule samplenames:
	output:
		outdir + "/logs/samples.names"
	message:
		"Creating list of sample names"
	run:
		with open(output[0], "w") as fout:
			for samplename in samplenames:
				_ = fout.write(samplename + "\n")	

rule bam_list:
    input: 
        bam = expand(bam_dir + "/{sample}.bam", sample = samplenames),
        bai = expand(bam_dir + "/{sample}.bam.bai", sample = samplenames)
    output:
        outdir + "/logs/samples.files"
    message:
        "Creating list of alignment files"
    run:
        with open(output[0], "w") as fout:
            for bamfile in input.bam:
                _ = fout.write(bamfile + "\n")

rule call_variants:
    input:
        bam = expand(bam_dir + "/{sample}.bam", sample = samplenames),
        bai = expand(bam_dir + "/{sample}.bam.bai", sample = samplenames),
        ref     = f"Assembly/{bn}",
        samples = outdir + "/logs/samples.files"
    output:
        temp(outdir + "/regions/{part}.vcf")
    message:
        "Calling variants: {wildcards.part}"
    log:
        outdir + "/logs/{part}.log"
    threads:
        2
    params:
        region = lambda wc: "-r " + regions[wc.part],
        ploidy = f"-p {ploidy}",
        populations = '' if groupings is None else f"--populations {groupings}",
        extra = extra
    shell:
        "freebayes -f {input.ref} -L {input.samples} {params} | bcftools sort - --output {output} 2> {log}"

rule vcf_list:
    output:
        outdir + "/logs/vcf.files"
    message:
        "Creating list of region-specific vcf files"
    run:
        with open(output[0], "w") as fout:
            for vcf in _regions:
                _ = fout.write(f"{outdir}/regions/{vcf}.vcf" + "\n")   

rule merge_vcfs:
    input:
        vcfs = expand(outdir + "/regions/{part}.vcf", part = _regions),
        filelist = outdir + "/logs/vcf.files"
    output:
        bcf = outdir + "/variants.raw.bcf",
        idx = outdir + "/variants.raw.bcf.csi"
    message:
        "Combining vcfs into a single file"
    log:
        "logs/concat.log"
    threads:
        50
    shell:  
        "bcftools concat -f {input.filelist} --threads {threads} --rm-dups -Ob --naive --write-index > {output.bcf} 2> {log}"

rule normalize_bcf:
	input: 
		genome  = f"Assembly/{bn}",
		bcf     = outdir + "/variants.raw.bcf"
	output:
		bcf     = outdir + "/variants.normalized.bcf",
		idx     = outdir + "/variants.normalized.bcf.csi",
	message: 
		"Normalizing the called variants"
	threads: 2
	shell:
		"""
		bcftools norm -d none -f {input.genome} {input.bcf} | bcftools norm -m -any -N -Ob --write-index > {output.bcf}
		"""

rule variants_stats:
	input:
        genome  = f"Assembly/{bn}",
		bcf     = outdir + "/variants.{type}.bcf",
		idx     = outdir + "/variants.{type}.bcf.csi",
        samples = outdir + "/logs/samples.names"
	output:
		outdir + "/stats/variants.{type}.stats",
	message:
		"Calculating variant stats: variants.{wildcards.type}.bcf"
	shell:
		"""
		bcftools stats -S {input.samples} --fasta-ref {input.genome} {input.bcf} > {output}
		"""

rule bcfreport:
    input:
        outdir + "/stats/variants.{type}.stats"
    output:
        outdir + "/stats/variants.{type}.html"
    message:
        "Generating bcftools report: variants.{wildcards.type}.bcf"
    script:
        "reportBcftools.Rmd"

rule log_runtime:
    output:
        outdir + "/logs/harpy.variants.log"
    message:
        "Creating record of relevant runtime parameters: {output}"
    params:
        ploidy = f"-p {ploidy}",
        populations = '' if groupings is None else f"--populations {groupings}",
        extra = extra
    run:
        with open([output[0]], "w") as f:
        	_ = f.write("The harpy variants module ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write(f"The directory with alignments: {bam_dir}\n")
            _ = f.write(f"Size of intervals to split genome for variant calling: {chunksize}\n")
            _ = f.write("The freebayes parameters:\n")
            _ = f.write(f"\tfreebayes -f GENOME -L samples.list -r REGION {params} | bcftools sort -\n")
            _ = f.write("The variants identified in the intervals were merged into the final variant file using:\n")
            _ = f.write("\tbcftools concat -f vcf.list --rm-dups --naive\t")
            _ = f.write("The variants were normalized using:\n")
			_ = f.write("\tbcftools norm -d none | bcftools norm -m -any -N -Ob")

rule all:
    input: 
        outdir + "/logs/harpy.variants.log",
        expand(outdir + "/variants.{file}.bcf",        file = ["raw", "normalized"]),
        expand(outdir + "/stats/variants.{file}.html", file = ["raw", "normalized"])
    message:
        "Variant calling is complete!"
    default_target: True