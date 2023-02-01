import os

bam_dir = config["seq_directory"]
genomefile = config["genomefile"]
groupings = config["groupings"]
ploidy = config["ploidy"]
samplenames = config["samplenames"]
extra = config["extra"]

def contignames(infile):
    with open(infile) as f:
        lines = [line.rstrip().split("\t")[0] for line in f]
    return lines

contigs = contignames(genomefile + ".fai")

rule index_alignments:
    input: bam_dir + "/{sample}.bam"
    output: bam_dir + "/{sample}.bam.bai"
    message: "Indexing barcodes: {wildcards.sample}"
    shell:
        """
        sambamba index {input} {output}
        """

rule split_contigs:
    input: f"{genomefile}.fai"
    output: temp(expand("Variants/mpileup/regions/{part}", part = contigs))
    message: "Separating {input} by contig for parallelization later"
    run:
        with open(input[0]) as f:
            cpath = "Variants/mpileup/regions"
            for line in f:
                contig = line.rstrip().split("\t")[0]
                with open(f"{cpath}/{contig}", "w") as fout:
                    gremlin = fout.write(f"{contig}\n")

rule bam_list:
    input: 
        bam = expand(bam_dir + "/{sample}.bam", sample = samplenames),
        bai = expand(bam_dir + "/{sample}.bam.bai", sample = samplenames)
    output: "Variants/mpileup/samples.list"
    message: "Creating list of alignment files"
    run:
        with open(output[0], "w") as fout:
            for bamfile in input.bam:
                fout.write(bamfile + "\n")

rule mpileup:
    input:
        bamlist = "Variants/mpileup/samples.list",
        genome = genomefile,
        region = "Variants/mpileup/regions/{part}"
    output: pipe("Variants/mpileup/{part}.mp.bcf")
    message: "Finding variants: {wildcards.part}"
    params:
        extra = extra
    shell:
        """
        bcftools mpileup --fasta-ref {input.genome} {params} --regions {wildcards.part} --bam-list {input.bamlist} --annotate AD --output-type b > {output} 2> /dev/null
        """

rule call_genotypes:
    input: "Variants/mpileup/{part}.mp.bcf"
    output: temp("Variants/mpileup/{part}.bcf")
    message: "Calling genotypes: {wildcards.part}"
    threads: 1
    params: 
        groupsamples = '' if groupings == 'none' else "--group-samples " + groupings,
        ploidy = f"--ploidy {ploidy}"
    shell:
        """
        bcftools call --multiallelic-caller {params} --variants-only --output-type b {input} | bcftools sort - --output {output} 2> /dev/null
        """

rule index_bcf:
    input: 
        bcf = "Variants/mpileup/{part}.bcf",
        samplelist = "Variants/mpileup/samples.list"
    output: temp("Variants/mpileup/{part}.bcf.csi")
    log: "Variants/mpileup/stats/{part}.stats"
    message: "Indexing: {wildcards.part}"
    threads: 4  
    shell:
        """
        bcftools index --threads {threads} --output {output} {input.bcf}
        bcftools stats {input} -S {input.samplelist} > {log}
        """

rule combine_bcfs:
    input: 
        bcf = expand("Variants/mpileup/{part}.bcf", part = contigs),
        idx = expand("Variants/mpileup/{part}.bcf.csi", part = contigs)
    output: 
        bcf = "Variants/mpileup/variants.raw.bcf",
        idx = "Variants/mpileup/variants.raw.bcf.csi",
        stats = "Variants/mpileup/variants.raw.stats"
    message: "Merging sample BCFs into: {output}"
    threads: 50
    shell:
        """
        bcftools concat --threads {threads} --output-type b --naive {input.bcf} > {output.bcf} 2> /dev/null
        bcftools index --output {output.idx} {output.bcf}
        bcftools stats {output.bcf} > {output.stats}
        """

rule bcfreport:
    input: "Variants/mpileup/variants.raw.stats"
    output: "Variants/mpileup/variants.raw.html"
    message: "Generating bcftools report: {output}"
    script: "../utilities/bcftoolsreport.Rmd"


rule all:
    input: 
        bcf = "Variants/mpileup/variants.raw.bcf",
        report = "Variants/mpileup/variants.raw.html"
    default_target: True