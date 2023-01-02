import os
import re

bam_dir = config["seq_directory"]
samplenames = config["samplenames"]
variantfile = config["variantfile"]

# Pull out the basename of the variant file
if variantfile.lower().endswith(".vcf"):
    ext = ".vcf"
elif variantfile.lower().endswith(".vcf.gz"):
    ext = ".vcf.gz"
elif variantfile.lower().endswith(".bcf"):
    ext = ".bcf"
else:
    print("ERROR: Supplied variant call file (" + variantfile + ") must end in one of [.vcf | .vcf.gz | .bcf]")
    exit(1)

# lazy method (in terms of effort) to remove the extension
variantbase = re.split(ext, os.path.basename(variantfile), flags = re.IGNORECASE)[0]

rule splitbysamplehet:
    input: 
        vcf = variantfile,
        bam = bam_dir + "/{sample}.bam"
    output: "Phasing/input/{sample}.het.bcf"
    message: "Extracting {wildcards.sample} from {input}"
    threads: 1
    shell:
        """
        bcftools view -s {wildcards.sample} -i 'INFO/INFO_SCORE >= 0.2' {input.vcf} |\\
        awk '/^#/;/CHROM/ {{OFS="\\t"}}; !/^#/ && $10~/^0\/1/' > {output}
        """

rule splitbysample:
    input: 
        vcf = variantfile,
        bam = bam_dir + "/{sample}.bam"
    output: "Phasing/input/{sample}.bcf"
    message: "Extracting {wildcards.sample} from {input}"
    threads: 1
    shell:
        """
        bcftools view -s {wildcards.sample} -i 'INFO/INFO_SCORE >= 0.2' {input.vcf} |\\
        awk '/^#/;/CHROM/ {{OFS="\\t"}}; !/^#/ &&  $10~/^0\/0/ {{$10="0|0:"substr($10,5);print $0}}; !/^#/ && $10~/^0\/1/; !/^#/ &&  $10~/^1\/1/ {{$10="1|1:"substr($10,5);print $0}}; !/^#/ {{print $0}}' > {output}
        """


#rule names:
#    input: expand("Phasing/input/{sample}.bcf", sample = samplenames)
#    output: 

rule extractHairs:
    input:
        vcf = "Phasing/input/{sample}.het.bcf",
        bam = bam_dir + "/{sample}.bam"
    output: "Phasing/extractHairs/{sample}.unlinked.frags"
    message: "Converting to compact fragment format: {wildcards.sample}"
    threads: 1
    shell:
        """
        extractHAIRS --10X 1 --nf 1 --bam {input.bam} --VCF {input.vcf} --out {output}
        """    

rule linkFragments:
    input: 
        bam = bam_dir + "/{sample}.bam",
        vcf = "Phasing/input/{sample}.het.bcf",
        fragments = "Phasing/extractHairs/{sample}.unlinked.frags"
    output: "Phasing/linkFragments/{sample}.linked.frags"
    message: "Linking fragments: {wildcards.sample}"
    params: d = 50000
    shell:
        """
        LinkFragments.py  --bam {input.bam} --VCF {input.vcf} --fragments {input.fragments} --out {output} -d {params};
        """


rule phaseBlocks:
    input:
        vcf = "Phasing/input/{sample}.het.bcf",
        fragments = "Phasing/linkFragments/{sample}.linked.frags"
    output: 
        blocks = "Phasing/phaseBlocks/{sample}.blocks",
        vcf = "Phasing/phaseBlocks/{sample}.blocks.phased.vcf"
    message: "Creating phased haplotype blocks: {wildcards.sample}"
    threads: 1
    shell:
        """
        HAPCUT2 --fragments {input.fragments} --vcf {input.vcf} --out {output.blocks} --nf 1 --threshold 30 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1
        """

rule createAnnotations:
    input: "Phasing/phaseBlocks/{sample}.blocks.phased.vcf"
    output: "Phasing/annotations/{sample}.annot.vcf.gz"
    message: "Creating annotation files: {wildcards.sample}"
    shell:
        """
        bcftools query -f "%CHROM\\t%POS[\\t%GT\\t%PS\\t%PQ\\t%PD]\\n" {input} | bgzip -c > {output}
        """

rule indexAnnotations:
    input: "Phasing/annotations/{sample}.annot.vcf.gz"
    output: "Phasing/annotations/{sample}.annot.vcf.gz.tbi"
    message: "Indexing {wildcards.sample}.annot.vcf.gz"
    shell:
        """
        tabix -b 2 -e 2 {input}     
        """

rule mergeAnnotations:
    input:
        annot = "Phasing/annotations/{sample}.annot.vcf.gz",
        orig = "Phasing/input/{sample}.bcf"
    output: 
        bcf = "Phasing/output/{sample}.phased.bcf",
        idx = "Phasing/output/{sample}.phased.bcf.csi"
    message: "Merging annotations: {wildcards.sample}"
    shell:
        """
        bcftools annotate -h add.hdr -a {input} -c CHROM,POS,FMT/GX,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT=1 |  awk '!/<ID=GX/' | sed 's/:GX:/:GT:/' | bcftools view - -o {output.bcf}; 
        bcftools index {output.bcf}
        """