from snakemake.utils import Paramspace
import pandas as pd
import subprocess

bam_dir = config["seq_directory"]
samplenames = config["samplenames"]
variantfile = config["variantfile"]

# declare a dataframe to be a paramspace
paramspace = Paramspace(pd.read_csv(config["paramfile"], sep="\t"), param_sep = "", filename_params="*")

def contignames(vcf):
    bcftools = subprocess.Popen(f"bcftools view --header-only {vcf}".split(), stdout = subprocess.PIPE)
    idcontigs = subprocess.Popen("grep contig=<ID".split(), stdin = bcftools.stdout, stdout = subprocess.PIPE)
    rm_suffix = subprocess.Popen("cut -d, -f1".split(), stdin = idcontigs.stdout, stdout = subprocess.PIPE)
    rm_prefix = subprocess.run(["sed", "s/##contig=<ID=//g"], stdin = rm_suffix.stdout, stdout = subprocess.PIPE)
    return sorted([chr for chr in rm_prefix.stdout.decode('utf-8').split()])

contigs = contignames(contigfile)
dict_cont = dict(zip(contigs, contigs))

rule bam_list:
    input: expand(bam_dir + "/{sample}.bam", sample = samplenames)
    output: "Imputation/input/samples.list"
    message: "Creating list of alignment files"
    benchmark: "Benchmark/Impute/filelist.txt"
    run:
        with open(output[0], "w") as fout:
            for bamfile in input:
                fout.write(f"{bamfile}\n")

rule prepare_biallelic_snps:
    input: variantfile
    output: pipe("Imputation/input/{part}.bisnp.bcf")
    message: "Keeping only biallelic SNPs from {wildcards.part}"
    params:
        lambda wc: dict_cont[wc.part]
    benchmark: "Benchmark/Impute/fileprep.{part}.txt"
    threads: 1
    shell:
        """
        bcftools view -m2 -M2 -v snps --regions {wildcards.part} --output-type b {input} > {output}
        """

#TODO investigate filter option
rule STITCH_format:
    input: "Imputation/input/{part}.bisnp.bcf"
    output: "Imputation/input/{part}.stitch"
    message: "Converting biallelic data to STITCH format: {wildcards.part}"
    benchmark: "Benchmark/Impute/stitchformat.{part}.txt"
    threads: 1
    params: 
        filters = "-i \'QUAL>20 && DP>10\'" if config["filtervcf"] else ""
    shell:
        """
        bcftools query {params} -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' {input} > {output}
        """

rule impute:
    input:
        bamlist = "Imputation/input/samples.list",
        infile = "Imputation/input/{part}.stitch"
    output:
        # format a wildcard pattern like "k{k}/s{s}/ngen{ngen}"
        # into a file path, with k, s, ngen being the columns of the data frame
        f"Imputation/{paramspace.wildcard_pattern}/contigs/" + "{part}/impute.vcf.gz"
    log: f"Imputation/{paramspace.wildcard_pattern}/contigs/" + "{part}/stitch.log"
    params:
        # automatically translate the wildcard values into an instance of the param space
        # in the form of a dict (here: {"k": ..., "s": ..., "ngen": ...})
        parameters = paramspace.instance
    message: "Running STITCH: {wildcards.part}\n  Parameters:\n  " + "{params.parameters}"
    benchmark: f"Benchmark/Impute/stitch.{paramspace.wildcard_pattern}" + ".{part}.txt"
    threads: 50
    script: "../utilities/stitch_impute.R"

rule samples_file:
    output: "Imputation/input/samples.names"
    message: "Creating file of sample names"
    threads: 1
    run:
        with open(output[0], "w") as fout:
            [fout.write(f"{i}\n") for i in samplenames]

rule index_vcf:
    input:
        vcf = "Imputation/{stitchparams}/contigs/{part}/impute.vcf.gz",
        samplelist = "Imputation/input/samples.names"
    output: 
        idx = "Imputation/{stitchparams}/contigs/{part}/impute.vcf.gz.tbi",
        stats = "Imputation/{stitchparams}/contigs/{part}/impute.stats"
    message: "Indexing: {wildcards.stitchparams}/{wildcards.part}"
    benchmark: "Benchmark/Impute/indexvcf.{stitchparams}.{part}.txt"
    threads: 1
    shell:
        """
        tabix {input.vcf}
        bcftools stats {input.vcf} -S {input.samplelist} > {output.stats}
        """

rule stitch_reports:
    input: "Imputation/{stitchparams}/contigs/{part}/impute.stats"
    output: "Imputation/{stitchparams}/contigs/{part}/{part}.report.html"
    message: "Generating STITCH report: {wildcards.part}"
    benchmark: "Benchmark/Impute/report.{stitchparams}.{part}.txt"
    threads: 1
    script: "../utilities/reportStitch.Rmd"

rule clean_stitch:
    input: "Imputation/{stitchparams}/contigs/{part}/{part}.report.html"
    message: "Cleaning up extra STITCH files"
    shell: 
        """
        rm -r Imputation/{wildcards.stitchparams}/contigs/{wildcards.part}/input
        rm -r Imputation/{wildcards.stitchparams}/contigs/{wildcards.part}/RData
        rm -r Imputation/{wildcards.stitchparams}/contigs/{wildcards.part}/plots
        """

rule merge_vcfs:
    input: 
        vcf = expand("Imputation/{{stitchparams}}/contigs/{part}/impute.vcf.gz", part = contigs),
        idx = expand("Imputation/{{stitchparams}}/contigs/{part}/impute.vcf.gz.tbi", part = contigs),
    output: "Imputation/{stitchparams}/variants.imputed.bcf"
    log: "Imputation/{stitchparams}/concat.log"
    message: "Merging VCFs: {wildcards.stitchparams}"
    benchmark: "Benchmark/Impute/mergevcf.{stitchparams}.txt"
    threads: 20
    shell:
        """
        bcftools concat --threads {threads} -o {output} --output-type b {input.vcf} 2> {log}
        """

rule stats:
    input:
        bcf = "Imputation/{stitchparams}/variants.imputed.bcf",
        samplelist = "Imputation/input/samples.names"
    output: "Imputation/{stitchparams}/variants.imputed.stats"
    message: "Indexing and calculating stats: {wildcards.stitchparams}/variants.imputed.bcf"
    benchmark: "Benchmark/Impute/mergestats.{stitchparams}.txt"
    shell:
        """
        bcftools index {input.bcf}
        bcftools stats {input.bcf} -S {input.samplelist} > {output}
        """

rule reports:
    input: "Imputation/{stitchparams}/variants.imputed.stats"
    output: "Imputation/{stitchparams}/variants.imputed.html"
    message: "Generating bcftools report: {output}"
    benchmark: "Benchmark/Impute/stitchreport.{stitchparams}.txt"
    script: "../utilities/reportBcftools.Rmd"

rule all:
    input: 
        bcf = expand("Imputation/{stitchparams}/variants.imputed.bcf", stitchparams=paramspace.instance_patterns),
        reports = expand("Imputation/{stitchparams}/variants.imputed.html", stitchparams=paramspace.instance_patterns),
        contigreports = expand("Imputation/{stitchparams}/contigs/{part}/{part}.report.html", stitchparams=paramspace.instance_patterns, part = contigs)
    default_target: True
    message: "Genotype imputation is complete!"