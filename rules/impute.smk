import os
import re
from snakemake.utils import Paramspace
import pandas as pd

# user specified configs
bam_dir = config["seq_directory"]
contigfile = config["contignames"]
#ncontigs = config["ncontigs"]
samplenames = config["samplenames"]
#model = config["method"]
#K = config["K"]
#S = config["S"]
#useBarcodes = str(config["useBarcodes"]).upper()
#nGenerations = config["nGenerations"]
variantfile = config["variantfile"]
#bx = "BX" if useBarcodes == "TRUE" else "noBX"
# declare a dataframe to be a paramspace
paramspace = Paramspace(pd.read_csv(config["paramfile"], sep="\t"), param_sep = "", filename_params="*")


# determine number of contigs from the contig file
def contigparts(contig_file):
    with open(contig_file, 'r') as fp:
        for ncontigs, line in enumerate(fp):
            pass
    ncontigs += 1
    return ncontigs

ncontigs = contigparts(contigfile)

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

# lazy method (in terms of effort) to remove the extension
variantbase = re.split(ext, os.path.basename(variantfile), flags = re.IGNORECASE)[0]

rule bam_list:
    input: expand(bam_dir + "/{sample}.bam", sample = samplenames)
    output: temp("Imputation/samples.list")
    message: "Creating list of alignment files"
    run:
        with open(output[0], "w") as fout:
            for bamfile in input:
                fout.write(bamfile + "\n")

rule split_contigs:
    input: contigfile
    output: expand("Imputation/contigs/contig.{part}", part = range(1, ncontigs + 1))
    message: "Splitting contig names for parallelization"
    shell:
        """
        awk '{{x="Imputation/contigs/contig."++i;}}{{print $1 > x;}}' {input}
        """

rule prepare_biallelic_snps:
    input: 
        vcf = variantfile,
        contig = "Imputation/contigs/contig.{part}"
    output: pipe("Imputation/input/" + variantbase + ".{part}.bisnp.bcf")
    message: "Keeping only biallelic SNPs from " + os.path.basename(variantfile) + ": contig {wildcards.part}"
    threads: 1
    shell:
        """
        bcftools view -m2 -M2 -v snps --regions $(cat {input.contig}) --output-type b {input.vcf} > {output}
        """

rule STITCH_format:
    input: "Imputation/input/" + variantbase + ".{part}.bisnp.bcf"
    output: "Imputation/input/" + variantbase + ".{part}"
    message: "Converting biallelic data to STITCH format: " + variantbase + ".{wildcards.part}"
    threads: 1
    shell:
        """
        bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {input} > {output}
        """

rule impute_search:
    input:
        bamlist = "Imputation/samples.list",
        infile = "Imputation/input/" + variantbase + ".{part}",
        chromosome = "Imputation/contigs/contig.{part}"
    output:
        # format a wildcard pattern like "k{k}/s{s}/ngen{ngen}"
        # into a file path, with k, s, ngen being the columns of the data frame
        f"Imputation/{paramspace.wildcard_pattern}/" + "contig{part}/contig{part}.impute.vcf.gz"
    log: f"Imputation/{paramspace.wildcard_pattern}/" + "contig{part}/contig{part}.stitch.log"
    params:
        # automatically translate the wildcard values into an instance of the param space
        # in the form of a dict (here: {"k": ..., "s": ..., "ngen": ...})
        parameters = paramspace.instance
    message: "Running STITCH: contig {wildcards.part}\n" + f"\n  parameters: {paramspace.instance}"
    threads: 50
    script: "../utilities/testparamspace.R"

#rule impute_genotypes:
#    input:
#        bamlist = "Imputation/samples.list",
#        infile = "Imputation/input/" + variantbase + ".{part}",
#        chromosome = "Imputation/contigs/contig.{part}"
#    output: "Imputation/" + model + "_K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "/contig{part}/K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "." + bx + model + ".vcf.gz"
#    log: "Imputation/" + model + "_K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "/contig{part}/K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "." + bx + model + ".log"
#    message: 
#        """
#        Running STITCH on contig {wildcards.part}
#        Parameters:
#            model: {params.model}
#            K: {params.K}
#            S: {params.S}
#            nGenerations: {params.nGenerations}
#            BX tags: {params.useBarcodes}
#        Execution progress logged to {log}
#        """
#    params:
#        model = model,
#        K = K,
#        S = S,
#        useBarcodes = useBarcodes,
#        nGenerations = nGenerations
#    threads: 50
#    script: "../utilities/stitch_impute.R"


rule testing:
    input: expand("Imputation/{params}/contig{part}/contig{part}.impute.vcf.gz", params=paramspace.instance_patterns, part = range(1, ncontigs + 1))
    default_target: True
    #input: expand(f"Imputation/{paramspace.wildcard_pattern}/" + "contig{part}/contig{part}.impute.vcf.gz", part = range(1, ncontigs + 1))

#rule testing:
#    input: expand("Imputation/" + model + "_K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "/contig{part}/K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "." + bx + model + ".vcf.gz", part = range(1, ncontigs + 1))
#    default_target: True


#rule vcf2bcf:
#    input: "Imputation/" + model + "_K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "/contig{part}.K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "." + bx + model + ".vcf"
#    output: "Imputation/" + model + "_K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "/contig{part}.K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "." + bx + model + ".bcf"
#    message: "Converting to BCF format: contig{part}" + ".K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "." + bx + "." + model
#    shell:
#        """
#        bcftools convert -Ob {input} | bcftools sort --output {output}
#        """
#
#rule index_bcf:
#    input: "Imputation/" + model + "_K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "/contig{part}.K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "." + bx + model + ".bcf"
#    output: temp("Imputation/" + model + "_K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "/contig{part}.K" + str(K) + "_S" + str(S) + "_nGen" + str(nGenerations) + "." + bx + model + ".bcf.csi")
#    message: "Indexing: {input}"
#    threads: 1
#    shell:
#        """
#        bcftools index --output {output} {input}
#        """
#
#rule merge_bcfs:
#    input: 
#        bcf = expand("VariantCall/leviathan/{sample}.bcf", sample = samplenames),
#        index = expand("VariantCall/leviathan/{sample}.bcf.csi", sample = samplenames)
#    output: "VariantCall/leviathan/variants.raw.bcf"
#    log: "VariantCall/leviathan/variants.raw.stats"
#    message: "Merging sample VCFs into single file: {output}"
#    #default_target: True
#    threads: 20
#    shell:
#        """
#        bcftools merge --threads {threads} -o {output} {input.bcf}
#        bcftools stats {output} > {log}
#        """
