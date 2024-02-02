from snakemake.utils import Paramspace
from rich import print as rprint
from rich.panel import Panel
import pandas as pd
import subprocess
import sys
import os

bam_dir     = config["seq_directory"]
samplenames = config["samplenames"]
variantfile = config["variantfile"]
paramfile   = config["paramfile"]
contigs     = config["contigs"]
skipreports = config["skipreports"]
# declare a dataframe to be the paramspace
paramspace  = Paramspace(pd.read_csv(paramfile, delim_whitespace = True).rename(columns=str.lower), param_sep = "", filename_params="*")

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

onerror:
    print("")
    rprint(
        Panel(
            f"The workflow has terminated due to an error. See the log file below for more details.",
            title = "[bold]harpy impute",
            title_align = "left",
            border_style = "red"
            ),
        file = sys.stderr
    )

onsuccess:
    print("")
    rprint(
        Panel(
            "The workflow has finished successfully! Find the results in [bold]Impute/[/bold]",
            title = "[bold]harpy impute",
            title_align = "left",
            border_style = "green"
            ),
        file = sys.stderr
    )

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
        "bcftools sort -Ob --write-index -o {output.bcf} {input} 2> {log}"

rule bam_list:
    input:
        expand(bam_dir + "/{sample}.bam", sample = samplenames)
    output:
        "Impute/input/samples.list"
    message:
        "Creating list of alignment files"
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
    threads: 
        3
    benchmark:
        ".Benchmark/Impute/fileprep.{part}.txt"
    message:
        "Converting data to biallelic STITCH format: {wildcards.part}"
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
        parameters = paramspace.instance,
        extra = ""
    conda:
        os.getcwd() + "/harpyenvs/r-env.yaml"
    benchmark:
        f".Benchmark/Impute/stitch.{paramspace.wildcard_pattern}" + ".{part}.txt"
    threads:
        50
    message: 
        "Performing imputation: {wildcards.part}\nmodel: {wildcards.model}\nuseBX: {wildcards.usebx}    \nbxLimit: {wildcards.bxlimit}\n    k: {wildcards.k}\n    s: {wildcards.s}\n nGen: {wildcards.ngen}"
    script:
        "stitch_impute.R"

rule index_vcf:
    input:
        vcf   = "Impute/{stitchparams}/contigs/{part}/{part}.vcf.gz"
    output: 
        idx   = "Impute/{stitchparams}/contigs/{part}/{part}.vcf.gz.tbi",
        stats = "Impute/{stitchparams}/contigs/{part}/{part}.stats"
    message:
        "Indexing: {wildcards.stitchparams}/{wildcards.part}"
    shell:
        """
        tabix {input.vcf}
        bcftools stats -s "-" {input.vcf} > {output.stats}
        """

rule collate_stitch_reports:
    input:
        "Impute/{stitchparams}/contigs/{part}/{part}.stats"
    output:
        "Impute/{stitchparams}/contigs/{part}/{part}.STITCH.html"
    message:
        "Generating STITCH report: {wildcards.part}"
    conda:
        os.getcwd() + "/harpyenvs/r-env.yaml"
    script:
        "report/StitchCollate.Rmd"

rule clean_stitch:
    input:
        "Impute/{stitchparams}/contigs/{part}/{part}.STITCH.html"
    output:
        temp("Impute/{stitchparams}/contigs/{part}/.cleaned")
    params:
        lambda wc: wc.get("stitchparams")
    priority:
        2
    message:
        "Cleaning up {wildcards.stitchparams}: {wildcards.part}"
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
        "Impute/{stitchparams}/variants.imputed.bcf"
    threads:
        workflow.cores
    message:
        "Merging VCFs: {wildcards.stitchparams}"
    shell:
        "bcftools concat --threads {threads} -O b -o {output} -f {input.files} 2> /dev/null"

rule index_merged:
    input:
        "Impute/{stitchparams}/variants.imputed.bcf"
    output:
        "Impute/{stitchparams}/variants.imputed.bcf.csi"
    message:
        "Indexing resulting file: {output}"
    shell:
        "bcftools index {input}"

rule stats:
    input:
        bcf = "Impute/{stitchparams}/variants.imputed.bcf",
        idx = "Impute/{stitchparams}/variants.imputed.bcf.csi"
    output:
        "Impute/{stitchparams}/reports/variants.imputed.stats"
    message:
        "Calculating stats: {wildcards.stitchparams}/variants.imputed.bcf"
    shell:
        """
        bcftools stats -s "-" {input.bcf} > {output}
        """

rule comparestats:
    input:
        orig    = "Impute/input/input.sorted.bcf",
        origidx = "Impute/input/input.sorted.bcf.csi",
        impute  = "Impute/{stitchparams}/variants.imputed.bcf",
        idx     = "Impute/{stitchparams}/variants.imputed.bcf.csi"
    output:
        compare = "Impute/{stitchparams}/reports/impute.compare.stats",
        info_sc = temp("Impute/{stitchparams}/reports/impute.infoscore")
    message:
        "Computing post-imputation stats: {wildcards.stitchparams}"
    shell:
        """
        bcftools stats -s "-" {input.orig} {input.impute} | grep \"GCTs\" > {output.compare}
        bcftools query -f '%CHROM\\t%POS\\t%INFO/INFO_SCORE\\n' {input.impute} > {output.info_sc}
        """

rule imputation_results_reports:
    input: 
        "Impute/{stitchparams}/reports/impute.compare.stats",
        "Impute/{stitchparams}/reports/impute.infoscore"
    output:
        "Impute/{stitchparams}/variants.imputed.html"
    params:
        lambda wc: wc.get("stitchparams")
    conda:
        os.getcwd() + "/harpyenvs/r-env.yaml"
    message:
        "Generating imputation success report: {output}"
    script:
        "report/Impute.Rmd"

rule log_runtime:
    output:
        "Impute/workflow/impute.workflow.summary"
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
            _ = f.write("\nThe STITCH parameters were governed by the rows of the input parameter table:\n")
            with open(config["paramfile"], "r") as f1:
                for line in f1:
                    _ = f.write("    " + line)
            _ = f.write("\nWithin R, STITCH was invoked with the following parameters:\n")
            _ = f.write(
                "    STITCH(\n" +
                "        method               = model,\n" +
                "        posfile              = posfile,\n" +
                "        bamlist              = bamlist,\n" +
                "        nCores               = ncores,\n" +
                "        nGen                 = ngen,\n" +
                "        chr                  = chr,\n" +
                "        K                    = k,\n" +
                "        S                    = s,\n" +
                "        use_bx_tag           = usebX,\n" +
                "        bxTagUpperLimit      = bxlimit,\n" +
                "        niterations          = 40,\n" +
                "        switchModelIteration = 39,\n" +
                "        splitReadIterations  = NA,\n" +
                "        outputdir            = outdir,\n" +
                "        output_filename      = outfile\n)\n"
            )
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]))

results = list()
results.append(expand("Impute/{stitchparams}/variants.imputed.bcf", stitchparams=paramspace.instance_patterns))
results.append(expand("Impute/{stitchparams}/contigs/{part}/{part}.STITCH.html", stitchparams=paramspace.instance_patterns, part = contigs))
results.append("Impute/workflow/impute.workflow.summary")

if not skipreports:
    results.append(expand("Impute/{stitchparams}/variants.imputed.html", stitchparams=paramspace.instance_patterns))

rule all:
    default_target: True
    input: 
        results
    message: 
        "Checking for expected workflow output"