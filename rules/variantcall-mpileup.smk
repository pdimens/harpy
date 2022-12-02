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
        samtools index {input}
        """

rule split_contigs:
    shell:
        """
        cat $fastafile | awk '$0 ~ "^>" {name=substr($0, 2); printf name"\t1\t"} $0 !~ "^>" {printf length($0)"\t"name"\n"}'
        """
rule mpileup:
    input:
        bam = bam_dir + "/{sample}" + ".bam",
        bai = bam_dir + "/{sample}" + ".bam.bai",
        genome = genomefile
    output: temp("VariantCall/{sample}.vcf")
    log:  "VariantCall/logs/{sample}.leviathan.log"
    message: "Calling variants: {wildcards.sample}"
    threads: 50
    shell:
        """
        bcftools mpileup \
            --fasta-ref {genome} \
            --regions $2 \
	        --annotate AD \
            --output-type u \
            --bam-list $1 | \

        """

shell:
    """
    bcftools call --multiallelic-caller \
        --group-samples POPMAPFILE \
        --variants-only \
        --output-type u - > ShadHap_all_${2}.bcf
    """