# user specified configs
bam_dir = config["seq_directory"]
genomefile = config["genome_file"]
groupings = config["groupings"]
n_regions = config["n_regions"]
ploidy = config["ploidy"]
samplenames = config["samplenames"] 

rule combine_bcfs:
    input: 
        bcf = expand("VariantCall/mpileup/region.{part}.bcf", part = range(1, n_regions + 1)),
        idx = expand("VariantCall/mpileup/region.{part}.bcf.csi", part = range(1, n_regions + 1))
    output: "VariantCall/mpileup/variants.raw.bcf"
    log: "VariantCall/mpileup/variants.raw.stats"
    message: "Merging sample BCFs into: {output}"
    default_target: True
    threads: 20
    shell:
        """
        bcftools concat --threads {threads} --output-type b --naive {input.bcf} > {output}
        bcftools stats {output} > {log}
        """

rule index_alignments:
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
    output: expand("VariantCall/mpileup/regions/region.{part}", part = range(1, n_regions + 1))
    message: "Separating {input} into regions for parallelization later"
    shell:
        """
        awk '{{x="VariantCall/mpileup/regions/region."++i;FS="\t"}}{{print $1 FS "1" FS $2 > x;}}' {input}
        """

rule bam_list:
    input: 
        bam = expand(bam_dir + "/{sample}.bam", sample = samplenames),
        bai = expand(bam_dir + "/{sample}.bam.bai", sample = samplenames)
    output: "VariantCall/mpileup/samples.list"
    message: "Creating list of alignment files"
    run:
        with open(output[0], "w") as fout:
            for bamfile in input.bam:
                fout.write(bamfile + "\n")

rule mpileup:
    input:
        bamlist = "VariantCall/mpileup/samples.list",
        genome = genomefile,
        region = "VariantCall/mpileup/regions/region.{part}"
    output: pipe("VariantCall/mpileup/region.{part}.mp.bcf")
    message: "Finding variants: region.{wildcards.part}"
    wildcard_constraints:
        part = "[0-9]*"
    threads: 1
    shell:
        """
        bcftools mpileup --fasta-ref {input.genome} --regions-file {input.region} --bam-list {input.bamlist} --annotate AD --output-type b > {output} 2> /dev/null
        """


rule call_genotypes:
        input: 
            bcf = "VariantCall/mpileup/region.{part}.mp.bcf"
        output: "VariantCall/mpileup/region.{part}.bcf"
        message: "Calling genotypes: region.{wildcards.part}"
        wildcard_constraints:
            part = "[0-9]*"
        threads: 1
        params: 
            groupsamples = '' if popfile == 'none' else "--group-samples " + groupings,
            ploidy = f"--ploidy {ploidy}"
        shell:
            """
            bcftools call --multiallelic-caller {params} --variants-only --output-type b {input.bcf} | bcftools sort - --output {output} 2> /dev/null
            """

rule index_bcf:
    input: "VariantCall/mpileup/region.{part}.bcf"
    output: temp("VariantCall/mpileup/region.{part}.bcf.csi")
    message: "Indexing: region.{wildcards.part}"
    wildcard_constraints:
        part = "[0-9]*"
    threads: 1    
    shell:
        """
        bcftools index --output {output} {input}
        """