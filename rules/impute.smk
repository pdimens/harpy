import re

# user specified configs
bam_dir = config["seq_directory"]
samplenames = config["samplenames"]
n_regions = config["n_regions"]
model = config["model"]
K = config["K"]
S = config["S"]
useBarcodes = config["useBarcodes"]
nGenerations = config["nGenerations"]
#contigs = config["contig_file"]
variantfile = config["variant_file"]
# Pull out the basename of the variant file
if variantfile.lower().endswith(".vcf"):
    variantbase = "".join(re.split(".vcf", variantfile, flags=re.IGNORECASE))
elif variantfile.lower().endswith(".vcf.gz"):
    variantbase = "".join(re.split(".vcf.gz", variantfile, flags=re.IGNORECASE))
elif variantfile.lower().endswith(".bcf"):
    variantbase = "".join(re.split(".bcf", variantfile, flags=re.IGNORECASE))
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

rule get_contigs:
    input: variantfile + ".positions"
    output: "Imputation/contigs.list"
    #output: temp(expand("Imputation/regions/region.{part}", part = range(1, n_regions + 1)))
    message: "Separating {input} into regions for STITCH parallelization"
    shell:
        """
        cut -d"\\t" -f1 {input} | sort | uniq > {output}
        """

rule split_contigs:
    input: "Imputation/contigs.list"
    output: temp(directory("Imputation/regions"))
    message: "Splitting contig names for parallelization"
    shell:
        """
        awk '{{x="{output}/region."++i;}}{{print $1 > x;}}' {input}
        """

rule impute_genotypes:
    input:
        bamlist = ,
        chromosome = 
    output: ""
    message: 
        """
        Running STITCH on {input.chromosome}
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
