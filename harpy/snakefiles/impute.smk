containerized: "docker://pdimens/harpy:latest"
import os
import logging

onstart:
    logger.logger.addHandler(logging.FileHandler(config["snakemake_log"]))
onsuccess:
    os.remove(logger.logfile)
onerror:
    os.remove(logger.logfile)
wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+",
    paramset = "[^/]+",
    contig = "[^/]+"

bamlist       = config["inputs"]["alignments"]
bamdict       = dict(zip(bamlist, bamlist))
variantfile   = config["inputs"]["variantfile"]
paramfile     = config["inputs"]["paramfile"]
biallelic     = config["inputs"]["biallelic_contigs"]
outdir        = config["output_directory"]
envdir        = os.path.join(os.getcwd(), outdir, "workflow", "envs")
skip_reports  = config["reports"]["skip"]
stitch_params = config["stitch_parameters"]
with open(biallelic, "r") as f:
    contigs = [line.rstrip() for line in f]

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
        bailist = collect("{bam}.bai", bam = bamlist)
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
        outdir + "/workflow/input/stitch/{contig}.stitch"
    threads: 
        3
    container:
        None
    shell:
        """
        bcftools view --types snps -M2 --regions {wildcards.contig} {input.bcf} |
            bcftools query -i '(STRLEN(REF)==1) & (STRLEN(ALT[0])==1) & (REF!="N")' -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n' > {output}
        """

rule impute:
    input:
        bamlist = outdir + "/workflow/input/samples.list",
        infile  = outdir + "/workflow/input/stitch/{contig}.stitch",
        bam = bamlist,
        bailist = collect("{bam}.bai", bam = bamlist)
    output:
        temp(directory(outdir + "/{paramset}/contigs/{contig}/plots")),
        temp(directory(outdir + "/{paramset}/contigs/{contig}/tmp")),
        temp(directory(outdir + "/{paramset}/contigs/{contig}/RData")),
        temp(directory(outdir + "/{paramset}/contigs/{contig}/input")),
        vcf = temp(outdir + "/{paramset}/contigs/{contig}/{contig}.vcf.gz")
    log:
        logfile = outdir + "/{paramset}/logs/{contig}.stitch.log"
    params:
        model   = lambda wc: stitch_params[wc.paramset]["model"],
        usebx   = lambda wc: stitch_params[wc.paramset]["usebx"],
        bxlimit = lambda wc: stitch_params[wc.paramset]["bxlimit"],
        k       = lambda wc: stitch_params[wc.paramset]["k"],
        s       = lambda wc: stitch_params[wc.paramset]["s"],
        ngen    = lambda wc: stitch_params[wc.paramset]["ngen"],
        extra   = config.get("stitch_extra", "")
    threads:
        workflow.cores - 1
    conda:
        f"{envdir}/stitch.yaml"
    script:
        "scripts/stitch_impute.R"

rule index_vcf:
    input:
        vcf   = outdir + "/{paramset}/contigs/{contig}/{contig}.vcf.gz"
    output:
        vcf   = outdir + "/{paramset}/contigs/{contig}.vcf.gz",
        idx   = outdir + "/{paramset}/contigs/{contig}.vcf.gz.tbi",
        stats = outdir + "/{paramset}/reports/data/contigs/{contig}.stats"
    container:
        None
    shell:
        """
        cp {input} {output.vcf}
        tabix {output.vcf}
        bcftools stats -s "-" {input.vcf} > {output.stats}
        """

rule contig_report:
    input:
        statsfile = outdir + "/{paramset}/reports/data/contigs/{contig}.stats",
        plotdir = outdir + "/{paramset}/contigs/{contig}/plots"
    output:
        outdir + "/{paramset}/reports/{contig}.{paramset}.html"
    log:
        logfile = outdir + "/{paramset}/logs/reports/{contig}.stitch.log"
    params:
        model   = lambda wc: stitch_params[wc.paramset]["model"],
        usebx   = lambda wc: stitch_params[wc.paramset]["usebx"],
        bxlimit = lambda wc: stitch_params[wc.paramset]["bxlimit"],
        k       = lambda wc: stitch_params[wc.paramset]["k"],
        s       = lambda wc: stitch_params[wc.paramset]["s"],
        ngen    = lambda wc: stitch_params[wc.paramset]["ngen"],
        extra   = config.get("stitch_extra", "")
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/stitch_collate.Rmd"

rule concat_list:
    input:
        bcf = collect(outdir + "/{{paramset}}/contigs/{contig}.vcf.gz", contig = contigs)
    output:
        temp(outdir + "/{paramset}/bcf.files")
    run:
        with open(output[0], "w") as fout:
            _ = fout.write("\n".join(input.bcf))

rule merge_vcf:
    priority: 100
    input:
        files = outdir + "/{paramset}/bcf.files",
        idx   = collect(outdir + "/{{paramset}}/contigs/{contig}.vcf.gz.tbi", contig = contigs)
    output:
        outdir + "/{paramset}/{paramset}.bcf"
    threads:
        workflow.cores
    container:
        None
    shell:
        "bcftools concat --threads {threads} -O b -o {output} -f {input.files} 2> /dev/null"

rule index_merged:
    input:
        outdir + "/{paramset}/{paramset}.bcf"
    output:
        outdir + "/{paramset}/{paramset}.bcf.csi"
    container:
        None
    shell:
        "bcftools index {input}"

rule general_stats:
    input:
        bcf = outdir + "/{paramset}/{paramset}.bcf",
        idx = outdir + "/{paramset}/{paramset}.bcf.csi"
    output:
        outdir + "/{paramset}/reports/data/impute.stats"
    container:
        None
    shell:
        "bcftools stats -s \"-\" {input.bcf} > {output}"

rule compare_stats:
    input:
        orig    = outdir + "/workflow/input/vcf/input.sorted.bcf",
        origidx = outdir + "/workflow/input/vcf/input.sorted.bcf.csi",
        impute  = outdir + "/{paramset}/{paramset}.bcf",
        idx     = outdir + "/{paramset}/{paramset}.bcf.csi"
    output:
        compare = outdir + "/{paramset}/reports/data/impute.compare.stats",
        info_sc = temp(outdir + "/{paramset}/reports/data/impute.infoscore")
    container:
        None
    shell:
        """
        bcftools stats -s "-" {input.orig} {input.impute} | grep \"GCTs\" > {output.compare}
        bcftools query -f '%CHROM\\t%POS\\t%INFO/INFO_SCORE\\n' {input.impute} > {output.info_sc}
        """

rule impute_reports:
    input: 
        comparison = outdir + "/{paramset}/reports/data/impute.compare.stats",
        infoscore = outdir + "/{paramset}/reports/data/impute.infoscore"
    output:
        outdir + "/{paramset}/reports/{paramset}.html"
    log:
        logfile = outdir + "/{paramset}/logs/reports/imputestats.log"
    params:
        paramname = lambda wc: wc.get("paramset"),
        model   = lambda wc: stitch_params[wc.paramset]["model"],
        usebx   = lambda wc: stitch_params[wc.paramset]["usebx"],
        bxlimit = lambda wc: stitch_params[wc.paramset]["bxlimit"],
        k       = lambda wc: stitch_params[wc.paramset]["k"],
        s       = lambda wc: stitch_params[wc.paramset]["s"],
        ngen    = lambda wc: stitch_params[wc.paramset]["ngen"],
        extra   = config.get("stitch_extra", "")
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/impute.Rmd"

rule workflow_summary:
    default_target: True
    input: 
        vcf = collect(outdir + "/{paramset}/{paramset}.bcf", paramset = list(stitch_params.keys())),
        agg_report = collect(outdir + "/{paramset}/reports/{paramset}.html", paramset = stitch_params.keys()) if not skip_reports else [],
        contig_report = collect(outdir + "/{paramset}/reports/{contig}.{paramset}.html", paramset = stitch_params.keys(), contig = contigs) if not skip_reports else [],
    run:
        paramfiletext = "\t".join(open(paramfile, "r").readlines())
        summary = ["The harpy impute workflow ran using these parameters:"]
        summary.append(f"The provided variant file: {variantfile}")
        preproc = "Preprocessing was performed with:\n"
        preproc += "\tbcftools view -M2 -v snps --regions CONTIG INFILE |\n"
        preproc += """\tbcftools query -i '(STRLEN(REF)==1) & (STRLEN(ALT[0])==1) & (REF!="N")' -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n'"""
        summary.append(preproc)
        stitchparam = f"The STITCH parameter file: {paramfile}\n"
        stitchparam += f"\t{paramfiletext}"
        summary.append(stitchparam)
        stitch = "Within R, STITCH was invoked with the following parameters:\n"
        stitch += "\tSTITCH(\n"
        stitch += "\t\tmethod = model,\n"
        stitch += "\t\tposfile = posfile,\n"
        stitch += "\t\tbamlist = bamlist,\n"
        stitch += "\t\tnCores = ncores,\n"
        stitch += "\t\tnGen = ngen,\n"
        stitch += "\t\tchr = chr,\n"
        stitch += "\t\tK = k,\n"
        stitch += "\t\tS = s,\n"
        stitch += "\t\tuse_bx_tag = usebx,\n"
        stitch += "\t\tbxTagUpperLimit = bxlimit,\n"
        stitch += "\t\tniterations = 40,\n"
        stitch += "\t\tswitchModelIteration = 39,\n"
        stitch += "\t\tsplitReadIterations = NA,\n"
        stitch += "\t\toutputdir = outdir,\n"
        stitch += "\t\toutput_filename = outfile\n\t)"
        stitchextra = "Additional STITCH parameters provided (overrides existing values above):\n"
        stitchextra += "\t" + config.get("stitch_extra", "None")
        summary.append(stitchextra)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['workflow_call']}"
        summary.append(sm)
        with open(outdir + "/workflow/impute.summary", "w") as f:
            f.write("\n\n".join(summary))
