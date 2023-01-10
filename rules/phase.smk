import os
import re

bam_dir = config["seq_directory"]
samplenames = config["samplenames"]
variantfile = config["variantfile"]

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
        awk '/^#/;/CHROM/ {{OFS="\\t"}}; !/^#/ &&  $10~/^0\\/0/ {{$10="0|0:"substr($10,5);print $0}}; !/^#/ && $10~/^0\\/1/; !/^#/ &&  $10~/^1\\/1/ {{$10="1|1:"substr($10,5);print $0}}; !/^#/ {{print $0}}' > {output}
        """

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
        LinkFragments.py  --bam {input.bam} --VCF {input.vcf} --fragments {input.fragments} --out {output} -d {params}
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
    output: "Phasing/annotations/{sample}.annot.bcf"
    message: "Creating annotation files: {wildcards.sample}"
    shell:
        """
        bcftools query -f "%CHROM\\t%POS[\\t%GT\\t%PS\\t%PQ\\t%PD]\\n" --output-type b {input} > {output}
        """

rule indexAnnotations:
    input: "Phasing/annotations/{sample}.annot.bcf"
    output: "Phasing/annotations/{sample}.annot.bcf.csi"
    message: "Indexing {wildcards.sample}.annot.bcf"
    shell:
        """
        bcftools index {input}     
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
        bcftools annotate -h add.hdr -a {input} -c CHROM,POS,FMT/GX,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT=1 |  awk '!/<ID=GX/' | sed 's/:GX:/:GT:/' | bcftools view - -o {output.bcf}
        bcftools index {output.bcf}
        """

rule mergeSamples:
    input:
        bcf = expand("Phasing/output/{sample}.phased.{ext}", sample = samplenames, ext = ["bcf", "bcf.csi"])
    output: "Phasing/output/variants.phased.bcf"
    default_target: True
    message: "Combinging samples into a single BCF file"
    threads: 30
    shell: "bcftools merge --threads {threads} --output-type b {input} > {output}"
