import re

# user specified configs
bam_dir = config["seq_directory"]
contigfile = config["contignames"]
ncontigs = config["ncontigs"]
samplenames = config["samplenames"]
model = config["model"]
K = config["K"]
S = config["S"]
useBarcodes = config["useBarcodes"]
nGenerations = config["nGenerations"]
#contigs = config["contig_file"]
variantfile = config["variant_file"]
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

rule bam_list:
    input: 
        bam = expand(bam_dir + "/{sample}.bam", sample = samplenames)
    output: temp("Imputation/samples.list")
    message: "Creating list of alignment files"
    run:
        with open(output[0], "w") as fout:
            for bamfile in input.bam:
                fout.write(bamfile + "\n")

rule prepare_input:
    input: 
        bcf = variantfile,
        regions = temp(expand("Impuatation/regions/region.{part}", part = range(1, n_regions + 1)))
    output: pipe(variantbase + ".biallelic.harpy.bcf")
    message: "Keeping only biallelic SNPs from {input}: contig {wildcards.part}"
    threads: 1
    shell:
        """
        bcftools view -m2 -M2 -v snps --output-type b {input} > {output}
        """

rule prepare_input:
    input: variantbase + ".biallelic.harpy.bcf"
    output: variantfile + ".positions"
    message: "Converting biallelic data to STITCH format"
    threads: 1
    shell:
        """
        bcftools query -f '%CHROM\\t%POS\\t$REF\\t%ALT\\n' {input} > {output}
        """

rule split_contigs:
    input: contigfile
    output: temp(expand("Imputation/contigs/contig.{part}", part = range(1, ncontigs + 1)))
    message: "Splitting contig names for parallelization"
    shell:
        """
        awk '{{x="Imputation/contigs/contig."++i;}}{{print $1 > x;}}' {input}
        """

#TODO fstring might not jive
bx = "BX" if useBarcodes.lower == "true" else "noBX"
rule impute_genotypes:
    input:
        bamlist = "Imputation/samples.list",
        chromosome = "Imputation/conrigs/contig.{part}"
    output: "Imputation/" + model + "_K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "/contig{part}.K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "." + bx + model + ".vcf"
    message: 
        """
        Running STITCH on contig {wildcards.part}
        Parameters:
            model: {params.model}
            K: {params.K}
            S: {params.S}
            nGenerations: {params.nGerenations}
            BX tags: {params.useBarcodes}
        """
    params:
        model = model,
        K = K,
        S = S,
        useBarcodes = useBarcodes,
        nGenerations = nGenerations
    threads: 50
    script: "utilities/stitch_impute.R"

rule vcf2bcf:
    input: "Imputation/" + model + "_K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "/contig{part}.K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "." + bx + model + ".vcf"
    output: "Imputation/" + model + "_K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "/contig{part}.K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "." + bx + model + ".bcf"
    message: "Converting to BCF format: contig{part}" + ".K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "." + bx + "." + model
    shell:
        """
        bcftools convert -Ob {input} | bcftools sort --output {output}
        """

rule index_bcf:
    input: "Imputation/" + model + "_K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "/contig{part}.K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "." + bx + model + ".bcf"
    output: temp("Imputation/" + model + "_K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "/contig{part}.K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "." + bx + model + ".bcf.csi")
    message: "Indexing: {input}"
    threads: 1
    shell:
        """
        bcftools index --output {output} {input}
        """

rule merge_bcfs:
    input: 
        bcf = expand("VariantCall/leviathan/{sample}.bcf", sample = samplenames),
        index = expand("VariantCall/leviathan/{sample}.bcf.csi", sample = samplenames)
    output: "VariantCall/leviathan/variants.raw.bcf"
    log: "VariantCall/leviathan/variants.raw.stats"
    message: "Merging sample VCFs into single file: {output}"
    default_target: True
    threads: 20
    shell:
        """
        bcftools merge --threads {threads} -o {output} {input.bcf}
        bcftools stats {output} > {log}
        """
