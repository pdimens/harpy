containerized: "docker://pdimens/harpy:latest"

import os
import pandas as pd
import logging
from snakemake.utils import Paramspace

onstart:
    logger.logger.addHandler(logging.FileHandler(config["snakemake_log"]))
onsuccess:
    os.remove(logger.logfile)
onerror:
    os.remove(logger.logfile)
wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

bamlist     = config["inputs"]["alignments"]
bamdict     = dict(zip(bamlist, bamlist))
variantfile = config["inputs"]["variantfile"]
paramfile   = config["inputs"]["paramfile"]
biallelic   = config["inputs"]["biallelic_contigs"]
outdir      = config["output_directory"]
envdir      = os.getcwd() + "/.harpy_envs"
skip_reports = config["skip_reports"]
paramspace  = Paramspace(pd.read_csv(paramfile, sep=r"\s+", skip_blank_lines=True).rename(columns=str.lower), param_sep = "", filename_params = ["k", "s", "ngen", "bxlimit"])
with open(biallelic, "r") as f_open:
    contigs = [i.rstrip() for i in f_open.readlines()]

rule sort_bcf:
    input:
        variantfile
    output:
        bcf = temp(f"{outdir}/workflow/input/vcf/input.sorted.bcf"),
        idx = temp(f"{outdir}/workflow/input/vcf/input.sorted.bcf.csi")
    log:
        f"{outdir}/logs/input.sort.log"
    container:
        None
    shell:
        "bcftools sort -Ob --write-index -o {output.bcf} {input} 2> {log}"

rule index_alignments:
    input:
        lambda wc: bamdict[wc.bam]
    output:
        "{bam}.bai"
    container:
        None
    shell:
        "samtools index {input}"

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
        temp(f"{outdir}/{paramspace.wildcard_pattern}/contigs/" + "{part}/{part}.vcf.gz"),
        temp(directory(f"{outdir}/{paramspace.wildcard_pattern}/contigs/" + "{part}/plots")),
        temp(directory(f"{outdir}/{paramspace.wildcard_pattern}/contigs/" + "{part}/tmp")),
        temp(directory(f"{outdir}/{paramspace.wildcard_pattern}/contigs/" + "{part}/RData")),
        temp(directory(f"{outdir}/{paramspace.wildcard_pattern}/contigs/" + "{part}/input"))
    log:
        f"{outdir}/{paramspace.wildcard_pattern}/logs/" + "{part}.stitch.log"
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
        workflow.cores - 1
    script:
        "scripts/stitch_impute.R"

rule index_vcf:
    input:
        vcf   = outdir + "/{stitchparams}/contigs/{part}/{part}.vcf.gz"
    output:
        vcf   = outdir + "/{stitchparams}/contigs/{part}.vcf.gz",
        idx   = outdir + "/{stitchparams}/contigs/{part}.vcf.gz.tbi",
        stats = outdir + "/{stitchparams}/reports/data/contigs/{part}.stats"
    wildcard_constraints:
        part = "[^/]+"
    container:
        None
    shell:
        """
        cp {input} {output.vcf}
        tabix {output.vcf}
        bcftools stats -s "-" {input.vcf} > {output.stats}
        """

rule stitch_reports:
    input:
        outdir + "/{stitchparams}/reports/data/contigs/{part}.stats",
        outdir + "/{stitchparams}/contigs/{part}/plots"
    output:
        outdir + "/{stitchparams}/reports/{part}.stitch.html"
    log:
        logfile = outdir + "/{stitchparams}/logs/reports/{part}.stitch.log"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/stitch_collate.Rmd"

rule concat_list:
    input:
        bcf = collect(outdir + "/{{stitchparams}}/contigs/{part}.vcf.gz", part = contigs)
    output:
        temp(outdir + "/{stitchparams}/bcf.files")
    run:
        with open(output[0], "w") as fout:
            _ = fout.write("\n".join(input.bcf))

rule merge_vcf:
    input:
        files = outdir + "/{stitchparams}/bcf.files",
        idx   = collect(outdir + "/{{stitchparams}}/contigs/{part}.vcf.gz.tbi", part = contigs)
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
        outdir + "/{stitchparams}/reports/data/impute.stats"
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
    log:
        logfile = outdir + "/{stitchparams}/logs/reports/imputestats.log"
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
        agg_report = collect(outdir + "/{stitchparams}/reports/variants.imputed.html", stitchparams=paramspace.instance_patterns) if not skip_reports else [],
        contig_report = collect(outdir + "/{stitchparams}/reports/{part}.stitch.html", stitchparams=paramspace.instance_patterns, part = contigs) if not skip_reports else [],
    run:
        paramfiletext = "\t".join(open(paramfile, "r").readlines())
        summary_template = f"""
The harpy impute workflow ran using these parameters:

The provided variant file: {variantfile}

Preprocessing was performed with:
    bcftools view -M2 -v snps --regions CONTIG INFILE |
    bcftools query -i '(STRLEN(REF)==1) & (STRLEN(ALT[0])==1) & (REF!="N")' -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n'

The STITCH parameter file: {paramfile}
    {paramfiletext}

Within R, STITCH was invoked with the following parameters:
    STITCH(
        method               = model,
        posfile              = posfile,
        bamlist              = bamlist,
        nCores               = ncores,
        nGen                 = ngen,
        chr                  = chr,
        K                    = k,
        S                    = s,
        use_bx_tag           = usebX,
        bxTagUpperLimit      = bxlimit,
        niterations          = 40,
        switchModelIteration = 39,
        splitReadIterations  = NA,
        outputdir            = outdir,
        output_filename      = outfile
    )

Additional STITCH parameters provided (overrides existing values above):\n")
    {config.get("extra", "None provided")}

The Snakemake workflow was called via command line:\n")
    {config["workflow_call"]}
"""
        with open(outdir + "/workflow/impute.summary", "w") as f:
            f.write(summary_template)
