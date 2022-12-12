# user specified configs
bam_dir = config["seq_directory"]
genomefile = config["genome_file"]
popfile = config["popfile"]

# Find the number of contigs in the genome fasta
n_regions = 0
with open(genomefile, "r") as fopen:
    while True:
        # Get next line from file
        line = fopen.readline(16)
        if line.startswith('>'):
            n_regions += 1
        # end of file is reached
        if not line:
            break

# Received from the harpy wrapper
samplenames = config["samplenames"] 

rule merge_vcfs:
    input: expand("VariantCall/{sample}.vcf", sample = samplenames)
    output: "VariantCall/variants.raw.bcf"
    log: "VariantCall/variants.raw.stats"
    message: "Merging sample BCFs into single file: {output}"
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

rule split_contigs:
    input: genomefile + ".fai"
    output: expand("VariantCall/regions/region.{part}", part = range(1, n_regions + 1))
    message: "Separating {input} into regions for parallelization later"
    params: n_regions
    shell:
        """
        awk '{{x="VariantCall/regions/region."++i;FS="\t"}}{{print $1 FS "0" FS $2 > x;}}' {input}
        """

rule bam_list:
    input: 
        bam = expand(bam_dir + "/{sample}.bam", sample = samplenames),
        bai = expand(bam_dir + "/{sample}.bam.bai", sample = samplenames)
    output: "VariantCall/samples.list"
    message: "Creating list of alignment files"
    run:
        with open(output[0], "w") as fout:
            for bamfile in input.bam:
                fout.write(bamfile + "\n")

rule mpileup:
    input:
        bamlist = "VariantCall/samples.list"
        genome = genomefile,
        region = "VariantCall/regions/region.{part}"
    output: pipe("VariantCall/region.{part}.mpileup.bcf")
    message: "Finding variants: region.{wildcards.part}"
    wildcard_constraints:
        part = "[0-9]*"
    threads: 1
    shell:
        """
        bcftools mpileup --fasta-ref {input.genome} --regions-file {input.region} --bam-list {input.bamlist} --annotate AD --output-type u > {output}
        """

rule call:
    input: 
        bcf = "VariantCall/region.{part}.mpileup.bcf",
        popmap = popfile
    output: "VariantCall/region.{part}.bcf"
    message: "Calling genotypes: region.{wildcard.part}"
    wildcard_constraints:
        part = "[0-9]*"
    threads: 1
    shell:
        """
        bcftools call --multiallelic-caller --group-samples {input.popfile} --variants-only --output-type u {input.bcf} > {output} 
        """