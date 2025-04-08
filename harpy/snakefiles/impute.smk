containerized: "docker://pdimens/harpy:latest"
import os
import re
import logging

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake_log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+",
    paramset = r"[^/]+",
    contig = r"[^/]+"

bamlist       = config["inputs"]["alignments"]
bamdict       = dict(zip(bamlist, bamlist))
variantfile   = config["inputs"]["variantfile"]
paramfile     = config["inputs"]["paramfile"]
region       = config.get("region", None)
envdir        = os.path.join(os.getcwd(), "workflow", "envs")
skip_reports  = config["reports"]["skip"]
stitch_params = config["stitch_parameters"]
stitch_extra  = config.get("stitch_extra", "None")
if region:
    contigs,positions = region.split(":")
    startpos,endpos,buffer = [int(i) for i in positions.split("-")]
    # remove the buffer to make it an htslib-style region
    region = f"{contigs}:{startpos}-{endpos}"
    # make the contig a list to fit with the existing workflow design
    contigs = [contigs]
else:
    biallelic = config["inputs"]["biallelic_contigs"]
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
        temp("workflow/input/vcf/input.sorted.bcf.csi"),
        bcf = temp("workflow/input/vcf/input.sorted.bcf")
    log:
        "logs/input.sort.log"
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
        collect("{bam}.bai", bam = bamlist),
        bam = bamlist,
    output:
        "workflow/input/samples.list"
    run:
        with open(output[0], "w") as fout:
            _ = [fout.write(f"{bamfile}\n") for bamfile in input.bam]

rule stitch_conversion:
    input:
        bcf = "workflow/input/vcf/input.sorted.bcf",
        idx = "workflow/input/vcf/input.sorted.bcf.csi"
    output:
        "workflow/input/stitch/{contig}.stitch"
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
        bamlist = "workflow/input/samples.list",
        infile  = "workflow/input/stitch/{contig}.stitch"
    output:
        temp(directory("{paramset}/contigs/{contig}/plots")),
        temp(directory("{paramset}/contigs/{contig}/RData")),
        temp(directory("{paramset}/contigs/{contig}/input")),
        temp("{paramset}/contigs/{contig}/{contig}.vcf.gz"),
        tmp = temp(directory("{paramset}/contigs/{contig}/tmp"))
    log:
        logfile = "{paramset}/logs/{contig}.stitch.log"
    params:
        chrom   = lambda wc: "--chr=" + wc.contig,
        start   = lambda wc: f"--regionStart={startpos}" if region else "",
        end     = lambda wc: f"--regionEnd={endpos}" if region else "",
        buffer  = lambda wc: f"--buffer={buffer}" if region else "",
        model   = lambda wc: "--method=" + stitch_params[wc.paramset]['model'],
        k       = lambda wc: f"--K={stitch_params[wc.paramset]['k']}",
        s       = lambda wc: f"--S={stitch_params[wc.paramset]['s']}",
        ngen    = lambda wc: f"--nGen={stitch_params[wc.paramset]['ngen']}",
        outdir  = lambda wc: "--outputdir=" + os.path.join(os.getcwd(), wc.paramset, "contigs", wc.contig),
        outfile = lambda wc: "--output_filename=" + f"{wc.contig}.vcf.gz",
        usebx   = lambda wc: "--use_bx_tag=" + str(stitch_params[wc.paramset]['usebx']).upper(),
        bxlimit = lambda wc: f"--bxTagUpperLimit={stitch_params[wc.paramset]['bxlimit']}",
        tmpdir  = lambda wc: "--tempdir=" + os.path.join(os.getcwd(), wc.paramset, "contigs", wc.contig, "tmp"),
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
        "{paramset}/contigs/{contig}/{contig}.vcf.gz"
    output:
        "{paramset}/contigs/{contig}.vcf.gz.tbi",
        vcf   = "{paramset}/contigs/{contig}.vcf.gz",
        stats = "{paramset}/reports/data/contigs/{contig}.stats"
    container:
        None
    shell:
        """
        cp {input} {output.vcf}
        tabix {output.vcf}
        bcftools stats -s "-" {input} > {output.stats}
        """

rule report_config:
    input:
        yaml = "workflow/report/_quarto.yml",
        scss = "workflow/report/_harpy.scss"
    output:
        yaml = temp("{paramset}/reports/_quarto.yml"),
        scss = temp("{paramset}/reports/_harpy.scss")
    run:
        import shutil
        for i,o in zip(input,output):
            shutil.copy(i,o)

rule contig_report:
    input:
        "{paramset}/reports/_harpy.scss",
        "{paramset}/reports/_quarto.yml",
        statsfile = "{paramset}/reports/data/contigs/{contig}.stats",
        plotdir = "{paramset}/contigs/{contig}/plots",
        qmd = "workflow/report/stitch_collate.qmd"
    output:
        report = "{paramset}/reports/{contig}.{paramset}.html",
        qmd = temp("{paramset}/reports/{contig}.{paramset}.qmd")
    log:
        logfile = "{paramset}/logs/reports/{contig}.stitch.log"
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
        cp -f {input.qmd} {output.qmd}
        STATS=$(realpath {input.statsfile})
        PLOTDIR=$(realpath {input.plotdir})
        quarto render {output.qmd} --log {log} --quiet -P statsfile:$STATS -P plotdir:$PLOTDIR {params}
        """

rule concat_list:
    input:
        collect("{{paramset}}/contigs/{contig}.vcf.gz", contig = contigs)
    output:
        temp("{paramset}/bcf.files")
    run:
        with open(output[0], "w") as fout:
            _ = fout.write("\n".join(input[0]))

if len(contigs) == 1:
    rule bcf_conversion:
        input:
            collect("{{paramset}}/contigs/{contig}.vcf.gz.tbi", contig = contigs),
            vcf = collect("{{paramset}}/contigs/{contig}.vcf.gz", contig = contigs)
        output:
            "{paramset}/{paramset}.bcf.csi",
            bcf = "{paramset}/{paramset}.bcf"
        container:
            None
        shell:
            "bcftools view -Ob --write-index -o {output.bcf} {input.vcf} 2> /dev/null"

else:
    rule merge_vcf:
        priority: 100
        input:
            collect("{{paramset}}/contigs/{contig}.vcf.gz.tbi", contig = contigs),
            files = "{paramset}/bcf.files"
        output:
            "{paramset}/{paramset}.bcf"
        threads:
            workflow.cores
        container:
            None
        shell:
            "bcftools concat --threads {threads} -Ob -o {output} -f {input.files} 2> /dev/null"

    rule index_merged:
        input:
            "{paramset}/{paramset}.bcf"
        output:
            "{paramset}/{paramset}.bcf.csi"
        container:
            None
        shell:
            "bcftools index {input}"

rule general_stats:
    input:
        "{paramset}/{paramset}.bcf.csi",
        bcf = "{paramset}/{paramset}.bcf"
    output:
        "{paramset}/reports/data/impute.stats"
    container:
        None
    shell:
        "bcftools stats -s \"-\" {input.bcf} > {output}"

rule extract_region:
    input:
        "workflow/input/vcf/input.sorted.bcf.csi",
        orig    = "workflow/input/vcf/input.sorted.bcf"
    output:
        temp("workflow/input/vcf/region.bcf.csi"),
        bcf = temp("workflow/input/vcf/region.bcf")
    params:
        f"-r {region}" if region else ""
    container:
        None
    shell:
        "bcftools view -Ob --write-index {params} -o {output.bcf} {input.orig}"
        
rule compare_stats:
    input:
        orig    = "workflow/input/vcf/input.sorted.bcf" if not region else "workflow/input/vcf/region.bcf",
        origidx = "workflow/input/vcf/input.sorted.bcf.csi" if not region else "workflow/input/vcf/region.bcf.csi",
        impute  = "{paramset}/{paramset}.bcf",
        idx     = "{paramset}/{paramset}.bcf.csi"
    output:
        compare = "{paramset}/reports/data/impute.compare.stats",
        info_sc = temp("{paramset}/reports/data/impute.infoscore")
    container:
        None
    shell:
        """
        bcftools stats -s "-" {input.orig} {input.impute} | grep \"GCTs\" > {output.compare}
        bcftools query -f '%CHROM\\t%POS\\t%INFO/INFO_SCORE\\n' {input.impute} > {output.info_sc}
        """

rule impute_reports:
    input:
        "{paramset}/reports/_quarto.yml",
        "{paramset}/reports/_harpy.scss",
        comparison = "{paramset}/reports/data/impute.compare.stats",
        infoscore = "{paramset}/reports/data/impute.infoscore",
        qmd = "workflow/report/impute.qmd"
    output:
        "{paramset}/reports/{paramset}.summary.html",
        qmd = temp("{paramset}/reports/{paramset}.summary.qmd")
    log:
        "{paramset}/logs/reports/imputestats.log"
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
        cp -f {input.qmd} {output.qmd}
        COMPARE=$(realpath {input.comparison})
        INFOSCORE=$(realpath {input.infoscore})
        quarto render {output.qmd} --log {log} --quiet -P compare:$COMPARE -P info:$INFOSCORE {params}
        """

rule workflow_summary:
    default_target: True
    input: 
        vcf = collect("{paramset}/{paramset}.bcf", paramset = list(stitch_params.keys())),
        agg_report = collect("{paramset}/reports/{paramset}.summary.html", paramset = stitch_params.keys()) if not skip_reports else [],
        contig_report = collect("{paramset}/reports/{contig}.{paramset}.html", paramset = stitch_params.keys(), contig = contigs) if not skip_reports else [],
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
        stitch += f"\t\tregionStart = {startpos},\n" if region else ""
        stitch += f"\t\tregionEnd = {endpos},\n" if region else ""
        stitch += f"\t\tbuffer = {buffer},\n" if region else ""
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
        sm += f"\t{config['snakemake_command']}"
        summary.append(sm)
        with open("workflow/impute.summary", "w") as f:
            f.write("\n\n".join(summary))
