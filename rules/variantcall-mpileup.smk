import os

bam_dir = config["seq_directory"]
genomefile = config["genomefile"]
groupings = config["groupings"]
ploidy = config["ploidy"]
samplenames = config["samplenames"]

def contignames(infile):
    with open(infile) as f:
        lines = [line.rstrip().split("\t")[0] for line in f]
    return lines

if not os.path.exists(f"{genomefile}.fai"):
    bn = os.path.basename(genomefile)
    print(f"{bn}.fai not found, indexing {bn} with samtools faidx")
    subprocess.run(["samtools","faidx", genomefile])

contigs = contignames(genomefile + ".fai")


rule combine_bcfs:
    input: 
        bcf = expand("Variants/mpileup/{part}.bcf", part = contigs),
        idx = expand("Variants/mpileup/{part}.bcf.csi", part = contigs)
    output: 
        bcf = "Variants/mpileup/variants.raw.bcf",
        idx = "Variants/mpileup/variants.raw.bcf.csi"
    log: report("Variants/mpileup/variants.raw.stats")
    message: "Merging sample BCFs into: {output}"
    default_target: True
    threads: 50
    shell:
        """
        bcftools concat --threads {threads} --output-type b --naive {input.bcf} > {output.bcf} 2> /dev/null
        bcftools index --output {output.idx} {output.bcf}
        bcftools stats {output.bcf} > {log}
        """

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
    shell:
        """
        awk '{{print $1 > "Variants/mpileup/regions/"$1;}}' {input}
        """

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
    shell:
        """
        bcftools mpileup --fasta-ref {input.genome} --regions {wildcards.part} --bam-list {input.bamlist} --annotate AD --output-type b > {output} 2> /dev/null
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