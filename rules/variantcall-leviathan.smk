# user specified configs
bam_dir = config["seq_directory"]
genomefile = config["genome_file"]
# Received from the harpy wrapper
samplenames = config["samplenames"] 

rule merge_vcfs:
    input: expand("VariantCall/{sample}.vcf.gz", sample = samplenames)
    output: "VariantCall/variants.raw.bcf"
    log: "VariantCall/variants.raw.stats"
    message: "Merging sample VCFs into single file: {output}"
    default_target: True
    threads: 20
    shell:
        """
        bcftools merge --threads {threads} -o {output} {input}
        bcftools stats {output} > {log}
        """

rule index_alignment:
    input: bam_dir + "/{sample}.bam"
    output: bam_dir + "/{sample}.bam.bai"
    message: "Indexing barcodes: {wildcards.sample}"
    threads: 1
    shell:
        """
        sambamba index {input} {output}
        """

rule index_barcode:
    input: 
        bam = bam_dir + "/{sample}.bam",
        bai = bam_dir + "/{sample}.bam.bai"
    output: temp("VariantCall/{sample}.bci")
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
    output: temp("VariantCall/{sample}.vcf")
    log:  "VariantCall/logs/{sample}.leviathan.log"
    message: "Calling variants: {wildcards.sample}"
    threads: 10
    shell:
        """
        LEVIATHAN -t {threads} -b {input.bam} -i {input.bc_idx} -g {input.genome} -o {output} 2> {log}
        """

rule compress_vcf:
    input: "VariantCall/{sample}.vcf"
    output: temp("VariantCall/{sample}.vcf.gz")
    message: "Compressing: {input}"
    threads: 5
    shell:        
        """
        bgzip --threads {threads} --stdout --reindex {input} > {output}
        """