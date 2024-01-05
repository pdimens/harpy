from snakemake.utils import Paramspace
import pandas as pd
import subprocess
import sys
import os

bam_dir     = config["seq_directory"]
samplenames = config["samplenames"]
variantfile = config["variantfile"]
paramfile   = config["paramfile"]
contigs     = config["contigs"]
# declare a dataframe to be the paramspace
paramspace  = Paramspace(pd.read_csv(paramfile, delim_whitespace = True).rename(columns=str.lower), param_sep = "", filename_params="*")

rule sort_bcf:
    input:
        variantfile
    output:
        bcf = temp("Impute/input/input.sorted.bcf"),
        idx = temp("Impute/input/input.sorted.bcf.csi")
    log:
        "Impute/input/input.sorted.log"
    message:
        "Sorting input variant call file"
    shell:
        """
        #bcftools sort -Ob --write-index {input} > {output.bcf} 2> {log}
        bcftools sort -Ob {input} > {output.bcf} 2> {log}
        bcftools index {output.bcf}
        """

rule bam_list:
    input:
        expand(bam_dir + "/{sample}.bam", sample = samplenames)
    output:
        "Impute/input/samples.list"
    message:
        "Creating list of alignment files"
    benchmark:
        ".Benchmark/Impute/filelist.txt"
    run:
        with open(output[0], "w") as fout:
            _ = [fout.write(f"{bamfile}\n") for bamfile in input]

rule samples_file:
    output:
        "Impute/input/samples.names"
    message:
        "Creating file of sample names"
    run:
        with open(output[0], "w") as fout:
            _ = [fout.write(f"{i}\n") for i in samplenames]

rule convert2stitch:
    input:
        "Impute/input/input.sorted.bcf"
    output:
        "Impute/input/{part}.stitch"
    message:
        "Converting data to biallelic STITCH format: {wildcards.part}"
    #params:
        #filters = "-i \'QUAL>20 && DP>10\'" if config["filtervcf"] else ""
    benchmark:
        ".Benchmark/Impute/fileprep.{part}.txt"
    threads: 3
    shell:
        """
        bcftools view --types snps -M2 --regions {wildcards.part} {input} |
        bcftools query -i '(STRLEN(REF)==1) & (STRLEN(ALT[0])==1) & (REF!="N")' -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' > {output}
        """

rule impute:
    input:
        bamlist = "Impute/input/samples.list",
        infile  = "Impute/input/{part}.stitch"
    output:
        # format a wildcard pattern like "k{k}/s{s}/ngen{ngen}"
        # into a file path, with k, s, ngen being the columns of the data frame
        f"Impute/{paramspace.wildcard_pattern}/contigs/" + "{part}/{part}.vcf.gz"
    log:
        f"Impute/{paramspace.wildcard_pattern}/contigs/" + "{part}/{part}.log"
    params:
        # automatically translate the wildcard values into an instance of the param space
        # in the form of a dict (here: {"k": ..., "s": ..., "ngen": ...})
        parameters = paramspace.instance
    message: 
        "Performing imputation: {wildcards.part}\nmodel: {wildcards.model}\nuseBX: {wildcards.usebx}    \nbxLimit: {wildcards.bxlimit}\n    k: {wildcards.k}\n    s: {wildcards.s}\n nGen: {wildcards.ngen}"
    benchmark:
        f".Benchmark/Impute/stitch.{paramspace.wildcard_pattern}" + ".{part}.txt"
    threads:
        50
    script:
        "stitch_impute.R"

rule index_vcf:
    input:
        vcf     = "Impute/{stitchparams}/contigs/{part}/{part}.vcf.gz"
        #samples = "Impute/input/samples.names"
    output: 
        idx        = "Impute/{stitchparams}/contigs/{part}/{part}.vcf.gz.tbi",
        stats      = "Impute/{stitchparams}/contigs/{part}/{part}.stats"
    message:
        "Indexing: {wildcards.stitchparams}/{wildcards.part}"
    benchmark:
        ".Benchmark/Impute/indexvcf.{stitchparams}.{part}.txt"
    shell:
        """
        tabix {input.vcf}
        bcftools stats -s "-" {input.vcf} > {output.stats}
        """

rule stitch_reports:
    input:
        "Impute/{stitchparams}/contigs/{part}/{part}.stats"
    output:
        "Impute/{stitchparams}/contigs/{part}/{part}.STITCH.html"
    message:
        "Generating STITCH report: {wildcards.part}"
    benchmark:
        ".Benchmark/Impute/report.{stitchparams}.{part}.txt"
    script:
        "reportImputeStitch.Rmd"

rule clean_stitch:
    input:
        "Impute/{stitchparams}/contigs/{part}/{part}.STITCH.html"
    output:
        temp("Impute/{stitchparams}/contigs/{part}/.cleaned")
    message:
        "Cleaning up {wildcards.stitchparams}: {wildcards.part}"
    priority:
        2
    params:
        lambda wc: wc.get("stitchparams")
    shell: 
        """
        rm -rf Impute/{params}/contigs/{wildcards.part}/input
        rm -rf Impute/{params}/contigs/{wildcards.part}/RData
        rm -rf Impute/{params}/contigs/{wildcards.part}/plots
        touch {output}
        """

rule concat_list:
    input:
        bcf = expand("Impute/{{stitchparams}}/contigs/{part}/{part}.vcf.gz", part = contigs)
    output:
        temp("Impute/{stitchparams}/bcf.files")
    message:
        "Creating list vcf files for concatenation"
    run:
        with open(output[0], "w") as fout:
            _ = fout.write("\n".join(input.bcf))

rule merge_vcfs:
    input:
        files = "Impute/{stitchparams}/bcf.files",
        idx   = expand("Impute/{{stitchparams}}/contigs/{part}/{part}.vcf.gz.tbi", part = contigs),
        clean = expand("Impute/{{stitchparams}}/contigs/{part}/.cleaned", part = contigs)
    output:
        bcf = "Impute/{stitchparams}/variants.imputed.bcf",
        bai = "Impute/{stitchparams}/variants.imputed.bcf.csi"
    message:
        "Merging VCFs: {wildcards.stitchparams}"
    benchmark:
        ".Benchmark/Impute/mergevcf.{stitchparams}.txt"
    threads:
        50
    shell:
        """
        #bcftools concat --threads {threads} --write-index -o {output.bcf}##idx##{output.bai} -f {input.files} 2> /dev/null
        bcftools concat --threads {threads} -o {output.bcf} -Ob -f {input.files} 2> /dev/null
        bcftools index {output.bcf}
        """

rule stats:
    input:
        bcf     = "Impute/{stitchparams}/variants.imputed.bcf",
        idx     = "Impute/{stitchparams}/variants.imputed.bcf.csi"
        #samples = "Impute/input/samples.names"
    output:
        "Impute/{stitchparams}/stats/variants.imputed.stats"
    message:
        "Calculating stats: {wildcards.stitchparams}/variants.imputed.bcf"
    benchmark:
        ".Benchmark/Impute/mergestats.{stitchparams}.txt"
    shell:
        """bcftools stats -s "-" {input.bcf} > {output}"""

rule comparestats:
    input:
        orig    = "Impute/input/input.sorted.bcf",
        origidx = "Impute/input/input.sorted.bcf.csi",
        impute  = "Impute/{stitchparams}/variants.imputed.bcf",
        idx     = "Impute/{stitchparams}/variants.imputed.bcf.csi"
        #samples = "Impute/input/samples.names"
    output:
        compare = "Impute/{stitchparams}/stats/impute.compare.stats",
        info_sc = temp("Impute/{stitchparams}/stats/impute.infoscore")
    message:
        "Computing post-imputation stats: {wildcards.stitchparams}"
    benchmark:
        ".Benchmark/Impute/mergestats.{stitchparams}.txt"
    shell:
        """
        bcftools stats -s "-" {input.orig} {input.impute} | grep \"GCTs\" > {output.compare}
        bcftools query -f '%CHROM\\t%POS\\t%INFO/INFO_SCORE\\n' {input.impute} > {output.info_sc}
        """

rule reports:
    input: 
        "Impute/{stitchparams}/stats/impute.compare.stats",
        "Impute/{stitchparams}/stats/impute.infoscore"
    output:
        "Impute/{stitchparams}/variants.imputed.html"
    message:
        "Generating imputation success report: {output}"
    params:
        lambda wc: wc.get("stitchparams")
    benchmark:
        ".Benchmark/Impute/stitchreport.{stitchparams}.txt"
    script:
        "reportImpute.Rmd"

rule log_runtime:
    output:
        "Impute/logs/harpy.impute.log"
    message:
        "Creating record of relevant runtime parameters: {output}"
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy impute module ran using these parameters:\n\n")
            _ = f.write(f"The provided variant file: {variantfile}\n")
            _ = f.write(f"The directory with alignments: {bam_dir}\n")
            _ = f.write("Preprocessing was performed with:\n")
            _ = f.write("    bcftools view -M2 -v snps --regions CONTIG INFILE |\n")
            _ = f.write("""    bcftools query -i '(STRLEN(REF)==1) & (STRLEN(ALT[0])==1) & (REF!="N")' -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n'\n""")
            _ = f.write("## STITCH imputation ##\n")
            _ = f.write("The STITCH parameters were governed by the rows of the input parameter table:\n")
            with open(config["paramfile"], "r") as f1:
                for line in f1:
                    _ = f.write(line)
            _ = f.write("\nWithin R, STITCH was invoked with the following parameters:\n")
            _ = f.write(
                "STITCH(\n" +
                "    method               = model,\n" +
                "    posfile              = posfile,\n" +
                "    bamlist              = bamlist,\n" +
                "    nCores               = ncores,\n" +
                "    nGen                 = ngen,\n" +
                "    chr                  = chr,\n" +
                "    K                    = k,\n" +
                "    S                    = s,\n" +
                "    use_bx_tag           = usebX,\n" +
                "    bxTagUpperLimit      = bxlimit,\n" +
                "    niterations          = 40,\n" +
                "    switchModelIteration = 39,\n" +
                "    splitReadIterations  = NA,\n" +
                "    outputdir            = outdir,\n" +
                "    output_filename      = outfile\n)\n"
            )

rule all:
    default_target: True
    input: 
        bcf     = expand("Impute/{stitchparams}/variants.imputed.bcf", stitchparams=paramspace.instance_patterns),
        reports = expand("Impute/{stitchparams}/variants.imputed.html", stitchparams=paramspace.instance_patterns),
        contigs = expand("Impute/{stitchparams}/contigs/{part}/{part}.STITCH.html", stitchparams=paramspace.instance_patterns, part = contigs),
        runlog  = "Impute/logs/harpy.impute.log"
    message: 
        "Genotype imputation is complete!"