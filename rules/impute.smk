import os
import re
from snakemake.utils import Paramspace
import pandas as pd

# user specified configs
bam_dir = config["seq_directory"]
contigfile = config["contignames"]
samplenames = config["samplenames"]
variantfile = config["variantfile"]
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
    params: 
        filters = "-i'QUAL>20 && DP>10'" if config["filtervcf"] else ""
    shell:
        """
        bcftools query {params} -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {input} > {output}
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
    message: "Running STITCH: contig {wildcards.part}\n  Parameters:\n  " + "{params.parameters}"
    wildcard_constraints:
        part = "[0-9]*"
    threads: 50
    script: "../utilities/testparamspace.R"

rule index_vcf:
    input: "Imputation/{stitchparams}/contig{part}/contig{part}.impute.vcf.gz"
    output: "Imputation/{stitchparams}/contig{part}/contig{part}.impute.vcf.gz.tbi"
    message: "Indexing: {wildcards.stitchparams}/contig{wildcards.part}"
    threads: 1
    wildcard_constraints:
        part = "[0-9]*"
    shell:
        """
        tabix {input}
        """

rule merge_vcfs:
    input: 
        vcf = expand("Imputation/{{stitchparams}}/contig{part}/contig{part}.impute.vcf.gz", part = range(1, ncontigs + 1)),
        idx = expand("Imputation/{{stitchparams}}/contig{part}/contig{part}.impute.vcf.gz.tbi", part = range(1, ncontigs + 1))
    output: 
        bcf = "Imputation/{stitchparams}/variants.imputed.bcf"
    log: 
        stats = "Imputation/{stitchparams}/variants.imputed.stats",
        concats = "Imputation/{stitchparams}/concat.log"
    message: "Merging VCFs: {wildcards.stitchparams}"
    threads: 20
    shell:
        """
        bcftools concat --threads {threads} -o {output} --output-type b {input.vcf} 2> {log.concats}
        bcftools stats {output} > {log.stats}
        """

rule reports:
    input: expand("Imputation/{stitchparams}/variants.imputed.bcf", stitchparams=paramspace.instance_patterns)
    output: "Imputation/report.html"
    message: "Generating report: {output}"
    default_target: True
    shell:
        """
        multiqc Imputation/model*useBX*/*.stats --force --quiet --no-data-dir --filename {output} 2> /dev/null
        """
    

