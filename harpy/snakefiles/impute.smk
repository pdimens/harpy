containerized: "docker://pdimens/harpy:latest"

import os
import re
import sys
import subprocess
import pandas as pd
import multiprocessing
import logging as pylogging
from snakemake.utils import Paramspace

bamlist     = config["inputs"]["alignments"]
variantfile = config["inputs"]["variantfile"]
paramfile   = config["inputs"]["paramfile"]
biallelic   = config["inputs"]["biallelic_contigs"]
outdir      = config["output_directory"]
envdir      = os.getcwd() + "/.harpy_envs"
skipreports = config["skip_reports"]
snakemake_log = config["snakemake_log"]
paramspace  = Paramspace(pd.read_csv(paramfile, sep=r"\s+", skip_blank_lines=True).rename(columns=str.lower), param_sep = "", filename_params = ["k", "s", "ngen", "bxlimit"])
with open(biallelic, "r") as f_open:
    contigs = [i.rstrip() for i in f_open.readlines()]

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

def sam_index(infile):
    """Use Samtools to index an input file, adding .bai to the end of the name"""
    if not os.path.exists(f"{infile}.bai"):
        subprocess.run(f"samtools index {infile} {infile}.bai".split())

onstart:
    extra_logfile_handler = pylogging.FileHandler(snakemake_log)
    logger.logger.addHandler(extra_logfile_handler)

rule sort_bcf:
    input:
        variantfile
    output:
        bcf = temp(f"{outdir}/workflow/input/vcf/input.sorted.bcf"),
        idx = temp(f"{outdir}/workflow/input/vcf/input.sorted.bcf.csi")
    log:
        f"{outdir}/workflow/input/vcf/input.sorted.log"
    container:
        None
    shell:
        "bcftools sort -Ob --write-index -o {output.bcf} {input} 2> {log}"

# not the ideal way of doing this, but it works
rule index_alignments:
    input:
        bamlist
    output:
        [f"{i}.bai" for i in bamlist]
    threads:
        workflow.cores
    run:
        with multiprocessing.Pool(processes=threads) as pool:
            pool.map(sam_index, input)

rule alignment_list:
    input:
        bam = bamlist,
        bailist = [f"{i}.bai" for i in bamlist]
    output:
        outdir + "/workflow/input/samples.list"
    run:
        with open(output[0], "w") as fout:
            _ = [fout.write(f"{bamfile}\n") for bamfile in input["bam"]]

rule stitch_conversion:
    input:
        bcf = f"{outdir}/workflow/input/vcf/input.sorted.bcf",
        idx = f"{outdir}/workflow/input/vcf/input.sorted.bcf.csi"
    output:
        outdir + "/workflow/input/stitch/{part}.stitch"
    threads: 
        3
    container:
        None
    shell:
        """
        bcftools view --types snps -M2 --regions {wildcards.part} {input.bcf} |
            bcftools query -i '(STRLEN(REF)==1) & (STRLEN(ALT[0])==1) & (REF!="N")' -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' > {output}
        """

rule impute:
    input:
        bamlist = outdir + "/workflow/input/samples.list",
        infile  = outdir + "/workflow/input/stitch/{part}.stitch"
    output:
        # format a wildcard pattern like "k{k}/s{s}/ngen{ngen}"
        # into a file path, with k, s, ngen being the columns of the data frame
        f"{outdir}/{paramspace.wildcard_pattern}/contigs/" + "{part}/{part}.vcf.gz"
    log:
        f"{outdir}/{paramspace.wildcard_pattern}/contigs/" + "{part}/{part}.log"
    params:
        # automatically translate the wildcard values into an instance of the param space
        # in the form of a dict (here: {"k": ..., "s": ..., "ngen": ...})
        parameters = paramspace.instance,
        extra = config.get("extra", "")
    conda:
        f"{envdir}/stitch.yaml"
    benchmark:
        f".Benchmark/{outdir}/stitch.{paramspace.wildcard_pattern}" + ".{part}.txt"
    threads:
        workflow.cores
    script:
        "scripts/stitch_impute.R"

rule index_vcf:
    input:
        vcf   = outdir + "/{stitchparams}/contigs/{part}/{part}.vcf.gz"
    output: 
        idx   = outdir + "/{stitchparams}/contigs/{part}/{part}.vcf.gz.tbi",
        stats = outdir + "/{stitchparams}/contigs/{part}/{part}.stats"
    container:
        None
    shell:
        """
        tabix {input.vcf}
        bcftools stats -s "-" {input.vcf} > {output.stats}
        """

rule stitch_reports:
    input:
        outdir + "/{stitchparams}/contigs/{part}/{part}.stats"
    output:
        outdir + "/{stitchparams}/contigs/{part}/{part}.STITCH.html"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/stitch_collate.Rmd"

rule clean_plots:
    input:
        outdir + "/{stitchparams}/contigs/{part}/{part}.STITCH.html"
    output:
        temp(touch(outdir + "/{stitchparams}/contigs/{part}/.plots.cleaned"))
    params:
        stitch = lambda wc: wc.get("stitchparams"),
        outdir = outdir
    container:
        None
    priority:
        2
    shell: 
        "rm -rf {params.outdir}/{params.stitch}/contigs/{wildcards.part}/plots"

rule concat_list:
    input:
        bcf = collect(outdir + "/{{stitchparams}}/contigs/{part}/{part}.vcf.gz", part = contigs)
    output:
        temp(outdir + "/{stitchparams}/bcf.files")
    run:
        with open(output[0], "w") as fout:
            _ = fout.write("\n".join(input.bcf))

rule merge_vcf:
    input:
        files = outdir + "/{stitchparams}/bcf.files",
        idx   = collect(outdir + "/{{stitchparams}}/contigs/{part}/{part}.vcf.gz.tbi", part = contigs)
    output:
        outdir + "/{stitchparams}/variants.imputed.bcf"
    threads:
        workflow.cores
    container:
        None
    shell:
        "bcftools concat --threads {threads} -O b -o {output} -f {input.files} 2> /dev/null"

rule index_merged:
    input:
        outdir + "/{stitchparams}/variants.imputed.bcf"
    output:
        outdir + "/{stitchparams}/variants.imputed.bcf.csi"
    container:
        None
    shell:
        "bcftools index {input}"

rule general_stats:
    input:
        bcf = outdir + "/{stitchparams}/variants.imputed.bcf",
        idx = outdir + "/{stitchparams}/variants.imputed.bcf.csi"
    output:
        outdir + "/{stitchparams}/reports/variants.imputed.stats"
    container:
        None
    shell:
        """
        bcftools stats -s "-" {input.bcf} > {output}
        """

rule compare_stats:
    input:
        orig    = outdir + "/workflow/input/vcf/input.sorted.bcf",
        origidx = outdir + "/workflow/input/vcf/input.sorted.bcf.csi",
        impute  = outdir + "/{stitchparams}/variants.imputed.bcf",
        idx     = outdir + "/{stitchparams}/variants.imputed.bcf.csi"
    output:
        compare = outdir + "/{stitchparams}/reports/data/impute.compare.stats",
        info_sc = temp(outdir + "/{stitchparams}/reports/data/impute.infoscore")
    container:
        None
    shell:
        """
        bcftools stats -s "-" {input.orig} {input.impute} | grep \"GCTs\" > {output.compare}
        bcftools query -f '%CHROM\\t%POS\\t%INFO/INFO_SCORE\\n' {input.impute} > {output.info_sc}
        """

rule impute_reports:
    input: 
        outdir + "/{stitchparams}/reports/data/impute.compare.stats",
        outdir + "/{stitchparams}/reports/data/impute.infoscore"
    output:
        outdir + "/{stitchparams}/reports/variants.imputed.html"
    params:
        lambda wc: wc.get("stitchparams")
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/impute.Rmd"


rule workflow_summary:
    default_target: True
    input: 
        vcf = collect(outdir + "/{stitchparams}/variants.imputed.bcf", stitchparams=paramspace.instance_patterns),
        agg_report = collect(outdir + "/{stitchparams}/reports/variants.imputed.html", stitchparams=paramspace.instance_patterns) if not skipreports else [],
        contig_report = collect(outdir + "/{stitchparams}/contigs/{part}/{part}.STITCH.html", stitchparams=paramspace.instance_patterns, part = contigs) if not skipreports else [],
        cleanedplots = collect(outdir + "/{stitchparams}/contigs/{part}/.plots.cleaned", stitchparams=paramspace.instance_patterns, part = contigs) if not skipreports else []
    message:
        "Summarizing the workflow: {output}"
    run:
        with open(outdir + "/workflow/impute.summary", "w") as f:
            _ = f.write("The harpy impute workflow ran using these parameters:\n\n")
            _ = f.write(f"The provided variant file: {variantfile}\n")
            _ = f.write("Preprocessing was performed with:\n")
            _ = f.write("    bcftools view -M2 -v snps --regions CONTIG INFILE |\n")
            _ = f.write("""    bcftools query -i '(STRLEN(REF)==1) & (STRLEN(ALT[0])==1) & (REF!="N")' -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n'\n""")
            _ = f.write("\nThe STITCH parameters were governed by the rows of the input parameter table:\n")
            with open(config["inputs"]["paramfile"], "r") as f1:
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
            _ = f.write("Additional STITCH parameters provided (overrides existing values above):\n")
            _ = f.write("    " + config.get("extra", "None provided") + "\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")