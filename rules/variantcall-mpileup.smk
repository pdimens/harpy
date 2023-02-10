import os

bam_dir = config["seq_directory"]
genomefile = config["genomefile"]
groupings = config.get("groupings", None)
ploidy = config["ploidy"]
samplenames = config["samplenames"]
extra = config.get("extra", "") 

def contignames(infile):
    with open(infile) as f:
        lines = [line.rstrip().split("\t")[0] for line in f]
    return lines

contigs = contignames(genomefile + ".fai")

rule index_alignments:
    input: bam_dir + "/{sample}.bam"
    output: bam_dir + "/{sample}.bam.bai"
    message: "Indexing barcodes: {wildcards.sample}"
    benchmark: "Benchmark/Variants/mpileup/indexbam.{sample}.txt"
    shell:
        """
        sambamba index {input} {output}
        """

rule split_contigs:
    input: f"{genomefile}.fai"
    output: temp(expand("Variants/mpileup/regions/{part}", part = contigs))
    message: "Separating {input} by contig for parallelization later"
    benchmark: "Benchmark/Variants/mpileup/splitcontigs.txt"
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
    benchmark: "Benchmark/Variants/mpileup/bamlist.txt"
    run:
        with open(output[0], "w") as fout:
            for bamfile in input.bam:
                _ = fout.write(bamfile + "\n")

rule mpileup:
    input:
        bamlist = "Variants/mpileup/samples.list",
        genome = genomefile,
        region = "Variants/mpileup/regions/{part}"
    output: pipe("Variants/mpileup/{part}.mp.bcf")
    message: "Finding variants: {wildcards.part}"
    log: "Variants/mpileup/logs/{part}.mpileup.log"
    benchmark: "Benchmark/Variants/mpileup/mpileup.{part}.txt"
    params:
        extra = extra
    shell:
        """
        bcftools mpileup --fasta-ref {input.genome} {params} --regions {wildcards.part} --bam-list {input.bamlist} --annotate AD --output-type b > {output} 2> {log}
        """

rule call_genotypes:
    input: "Variants/mpileup/{part}.mp.bcf"
    output: temp("Variants/mpileup/{part}.bcf")
    message: "Calling genotypes: {wildcards.part}"
    benchmark: "Benchmark/Variants/mpileup/call.{part}.txt"
    log: "Variants/mpileup/logs/{part}.call.log"
    threads: 2
    params: 
        groupsamples = '' if groupings is not None else f"--group-samples {groupings}",
        ploidy = f"--ploidy {ploidy}"
    shell:
        """
        bcftools call --multiallelic-caller {params} --variants-only --output-type b {input} | bcftools sort - --output {output} 2> /dev/null
        """

rule index_bcf:
    input: 
        bcf = "Variants/mpileup/{part}.bcf",
        samplelist = "Variants/mpileup/samples.list",
        genome = genomefile
    output: temp("Variants/mpileup/{part}.bcf.csi")
    log: "Variants/mpileup/stats/{part}.stats"
    message: "Indexing: {wildcards.part}"
    benchmark: "Benchmark/Variants/mpileup/indexbcf.{part}.txt"
    threads: 4
    shell:
        """
        bcftools index --threads {threads} --output {output} {input.bcf}
        bcftools stats -S {input.samplelist} --fasta-ref {input.genome} {input.bcf} > {log}
        """

rule combine_bcfs:
    input: 
        bcf = expand("Variants/mpileup/{part}.bcf", part = contigs),
        idx = expand("Variants/mpileup/{part}.bcf.csi", part = contigs),
        genome = genomefile
    output: 
        bcf = "Variants/mpileup/variants.raw.bcf",
        idx = "Variants/mpileup/variants.raw.bcf.csi",
        stats = "Variants/mpileup/variants.raw.stats"
    message: "Merging all BCFs into: {output.bcf}"
    benchmark: "Benchmark/Variants/mpileup/merge.txt"
    threads: 50
    shell:
        """
        bcftools concat --threads {threads} --output-type b --naive {input.bcf} > {output.bcf} 2> /dev/null
        bcftools index --output {output.idx} {output.bcf}
        bcftools stats --fasta-ref {input.genome} {output.bcf} > {output.stats}
        """

rule bcfreport:
    input: "Variants/mpileup/variants.raw.stats"
    output: "Variants/mpileup/variants.raw.html"
    message: "Generating bcftools report: {output}"
    benchmark: "Benchmark/Variants/mpileup/reports.txt"
    script: "../utilities/bcftoolsreport.Rmd"


rule all:
    input: 
        bcf = "Variants/mpileup/variants.raw.bcf",
        report = "Variants/mpileup/variants.raw.html"
    default_target: True