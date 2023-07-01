import os
import sys

bam_dir 	= config["seq_directory"]
genomefile 	= config["genomefile"]
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
    params:
        region = lambda wc: "-r " + regions[wc.part],
        ploidy = f"-p {ploidy}",
        populations = '' if groupings is None else f"--populations {groupings}",
        extra = extra
    shell:
        "freebayes -f {input.ref} -L {input.samples} {params} > {output} 2> {log}"

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
        outdir + "/variants.raw.bcf"
    message:
        "Combining vcfs into a single file"
    log:
        "logs/concat.log"
    threads:
        3
    shell:  
        "bcftools concat -f {input.filelist} --rm-dups | bcftools view -Ob > {output} 2> {log}"

rule index_bcf:
    input: 
        bcf = outdir + "/variants.raw.bcf",
        samples = outdir + "/logs/samples.names",
        genome  = f"Assembly/{genomefile}"
    output:
        outdir + "/variants.raw.bcf.csi"
    message:
        "Indexing: {input.bcf}"
    threads:
        4
    shell:
        """
        bcftools index --threads {threads} --output {output} {input.bcf}
        """

rule bcf_stats:
    input: 
        bcf = outdir + "/variants.raw.bcf",
        csi = outdir + "/variants.raw.bcf.csi",
        samples = outdir + "/logs/samples.names",
        genome  = f"Assembly/{genomefile}"
    output:
        outdir + "/stats/variants.raw.stats"
    message:
        "Calculating stats"
    benchmark:
        "Benchmark/Variants/mpileup/indexbcf.{part}.txt"
    shell:
        """
        bcftools stats -S {input.samples} --fasta-ref {input.genome} {input.bcf} > {output}
        """

rule bcfreport:
    input:
        outdir + "/stats/variants.raw.stats"
    output:
        outdir + "/stats/variants.raw.html"
    message:
        "Generating bcftools report: variants.raw.bcf"
    script:
        "reportBcftools.Rmd"

rule all:
    input: 
        outdir + "/variants.raw.bcf",
        outdir + "/stats/variants.raw.html"
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