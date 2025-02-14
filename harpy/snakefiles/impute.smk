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
stitch_extra  = config.get("stitch_extra", "None")
with open(biallelic, "r") as f:
    contigs = [line.rstrip() for line in f]

# instantiate static STITCH arguments
extraparams = {
    "--niterations" : 40,
    "--switchModelIteration" : 39,
    "--splitReadIterations" : "NA",
}
# update and overwrite with extraparms, if any
if stitch_extra != "None":
    for i in stitch_extra.split():
        arg,val = i.split("=")
        extraparams[arg] = val

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
        bamlist,
        collect("{bam}.bai", bam = bamlist),
        bamlist = outdir + "/workflow/input/samples.list",
        infile  = outdir + "/workflow/input/stitch/{contig}.stitch"
    output:
        temp(directory(outdir + "/{paramset}/contigs/{contig}/plots")),
        temp(directory(outdir + "/{paramset}/contigs/{contig}/RData")),
        temp(directory(outdir + "/{paramset}/contigs/{contig}/input")),
        temp(outdir + "/{paramset}/contigs/{contig}/{contig}.vcf.gz"),
        tmp = temp(directory(outdir + "/{paramset}/contigs/{contig}/tmp"))
    log:
        logfile = outdir + "/{paramset}/logs/{contig}.stitch.log"
    params:
        chrom   = lambda wc: "--chr=" + wc.contig,
        model   = lambda wc: "--method=" + stitch_params[wc.paramset]['model'],
        k       = lambda wc: f"--K={stitch_params[wc.paramset]['k']}",
        s       = lambda wc: f"--S={stitch_params[wc.paramset]['s']}",
        ngen    = lambda wc: f"--nGen={stitch_params[wc.paramset]['ngen']}",
        outdir  = lambda wc: "--outputdir=" + os.path.join(outdir, wc.paramset, "contigs", wc.contig),
        outfile = lambda wc: "--output_filename=" + f"{wc.contig}.vcf.gz",
        usebx   = lambda wc: "--use_bx_tag=" + str(stitch_params[wc.paramset]['usebx']).upper(),
        bxlimit = lambda wc: f"--bxTagUpperLimit={stitch_params[wc.paramset]['bxlimit']}",
        tmpdir  = lambda wc: "--tempdir=" + os.path.join(outdir, wc.paramset, "contigs", wc.contig, "tmp"),
        extra   = " ".join([f"{i}={j}" for i,j in extraparams.items()])
    threads:
        workflow.cores - 1
    conda:
        f"{envdir}/stitch.yaml"
    shell:
        """
        mkdir -p {output.tmp}
        STITCH.R --nCores={threads} --bamlist={input.bamlist} --posfile={input.infile} {params} 2> {log}
        """

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

rule report_config:
    input:
        yaml = f"{outdir}/workflow/report/_quarto.yml",
        scss = f"{outdir}/workflow/report/_harpy.scss"
    output:
        yaml = temp(f"{outdir}/{{paramset}}/reports/_quarto.yml"),
        scss = temp(f"{outdir}/{{paramset}}/reports/_harpy.scss")
    run:
        import shutil
        for i,o in zip(input,output):
            shutil.copy(i,o)

rule contig_report:
    input:
        f"{outdir}/{{paramset}}/reports/_quarto.yml",
        f"{outdir}/{{paramset}}/reports/_harpy.scss",
        statsfile = outdir + "/{paramset}/reports/data/contigs/{contig}.stats",
        plotdir = outdir + "/{paramset}/contigs/{contig}/plots",
        qmd = f"{outdir}/workflow/report/stitch_collate.qmd"
    output:
        report = outdir + "/{paramset}/reports/{contig}.{paramset}.html",
        qmd = temp(outdir + "/{paramset}/reports/{contig}.{paramset}.qmd")
    log:
        logfile = outdir + "/{paramset}/logs/reports/{contig}.stitch.log"
    params:
        params  = lambda wc: f"-P id:{wc.paramset}-{wc.contig}",
        model   = lambda wc: f"-P model:{stitch_params[wc.paramset]['model']}",
        usebx   = lambda wc: f"-P usebx:{stitch_params[wc.paramset]['usebx']}",
        bxlimit = lambda wc: f"-P bxlimit:{stitch_params[wc.paramset]['bxlimit']}",
        k       = lambda wc: f"-P k:{stitch_params[wc.paramset]['k']}",
        s       = lambda wc: f"-P s:{stitch_params[wc.paramset]['s']}",
        ngen    = lambda wc: f"-P ngen:{stitch_params[wc.paramset]['ngen']}",
        extra   = f"-P extra:{stitch_extra}"
    conda:
        f"{envdir}/r.yaml"
    shell:
        """
        cp {input.qmd} {output.qmd}
        STATS=$(realpath {input.statsfile})
        PLOTDIR=$(realpath {input.plotdir})
        quarto render {output.qmd} --log {log} --quiet -P statsfile:$STATS -P plotdir:$PLOTDIR {params}
        """

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
        f"{outdir}/{{paramset}}/reports/_quarto.yml",
        f"{outdir}/{{paramset}}/reports/_harpy.scss",
        comparison = outdir + "/{paramset}/reports/data/impute.compare.stats",
        infoscore = outdir + "/{paramset}/reports/data/impute.infoscore",
        qmd = f"{outdir}/workflow/report/impute.qmd"
    output:
        report = outdir + "/{paramset}/reports/{paramset}.summary.html",
        qmd = temp(outdir + "/{paramset}/reports/{paramset}.summary.qmd")
    log:
        outdir + "/{paramset}/logs/reports/imputestats.log"
    params:
        param   = lambda wc: f"-P id:{wc.paramset}",
        model   = lambda wc: f"-P model:{stitch_params[wc.paramset]['model']}",
        usebx   = lambda wc: f"-P usebx:{stitch_params[wc.paramset]['usebx']}",
        bxlimit = lambda wc: f"-P bxlimit:{stitch_params[wc.paramset]['bxlimit']}",
        k       = lambda wc: f"-P k:{stitch_params[wc.paramset]['k']}",
        s       = lambda wc: f"-P s:{stitch_params[wc.paramset]['s']}",
        ngen    = lambda wc: f"-P ngen:{stitch_params[wc.paramset]['ngen']}",
        extra   = f"-P extra:{stitch_extra}"
    conda:
        f"{envdir}/r.yaml"
    shell:
        """
        cp {input.qmd} {output.qmd}
        COMPARE=$(realpath {input.comparison})
        INFOSCORE=$(realpath {input.infoscore})
        quarto render {output.qmd} --log {log} --quiet -P compare:$COMPARE -P info:$INFOSCORE {params}
        """

rule workflow_summary:
    default_target: True
    input: 
        vcf = collect(outdir + "/{paramset}/{paramset}.bcf", paramset = list(stitch_params.keys())),
        agg_report = collect(outdir + "/{paramset}/reports/{paramset}.summary.html", paramset = stitch_params.keys()) if not skip_reports else [],
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
