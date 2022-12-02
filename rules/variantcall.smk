import os

# user specified configs
bam_dir = config["seq_directory"]
genomefile = config["genome_file"]
# Received from the harpy wrapper
samplenames = config["samplenames"] 

rule merge_vcfs:
    input: expand("VariantCall/{sample}.vcf", sample = samplenames)
    output: "VariantCall/variants.raw.bcf"
    log: "VariantCall/variants.raw.stats"
    message: "Merging sample VCFs into single file: {output}"
    default_target: True
    threads: 20
    shell:
        """
        bcftools merge --threads --output-type b -o {output} {input}
        bcftools stats {output} > {log}
        """

rule barcode_index:
    input: 
        bam = bam_dir + "/{sample}.bam",
        bai = bam_dir + "/{sample}.bam.bai"
    output: "VariantCall/{sample}.bci"
    message: "Indexing barcodes: {wildcards.sample}"
    threads: 4
    shell:
        """
        LRez index bam -p -b {input} -o {output}
        """

rule leviathan_variantcall:
    input:
        bam = bam_dir + "/{sample}" + ".bam",
        bai = bam_dir + "/{sample}" + ".bam.bai",
        bc_idx = "VariantCall/{sample}.bci",
        genome = genomefile
    output: "VariantCall/{sample}.vcf"
    message: "Calling variants: {wildcards.sample}"
    threads: 50
    shell:
        """
        LEVIATHAN -t {threads} -b {input.bam} -i {input.bc_idx} -g {input.genome} -o {output}      
        """