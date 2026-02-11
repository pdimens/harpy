import os
import re

wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+",
    paramset = r"[^/]+",
    contig = r"[^/]+"

VERSION=4.0
bamlist       = config["Inputs"]["alignments"]
bamdict       = dict(zip(bamlist, bamlist))
variantfile   = config["Inputs"]["vcf"]
paramfile     = config["Inputs"]["parameters"]
skip_reports  = config["Workflow"]["reports"]["skip"]
region        = config["Parameters"].get("region", None)
stitch_params = config["Parameters"]["stitch"]
stitch_extra  = config["Parameters"].get("extra", "None")
grid_size     = config["Parameters"]["grid-size"]
if region:
    contigs,positions = region.split(":")
    startpos,endpos,buffer = [int(i) for i in positions.split("-")]
    # remove the buffer to make it an htslib-style region
    region = f"{contigs}:{startpos}-{endpos}"
    # make the contig a list to fit with the existing workflow design
    contigs = [contigs]
else:
    biallelic = config["Inputs"]["biallelic-contigs"]
    with open(biallelic, "r") as f:
        contigs = [line.rstrip() for line in f]

# instantiate static STITCH arguments
extraparams = {
     **({"--gridWindowSize" : grid_size} if grid_size > 1 else {}),
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
        temp("{paramset}/contigs/{contig}/plots/alphaMat.all.png"),
        temp("{paramset}/contigs/{contig}/plots/alphaMat.normalized.png"),
        temp("{paramset}/contigs/{contig}/plots/hapSum_log.png"),
        temp("{paramset}/contigs/{contig}/plots/hapSum.png"),
        temp("{paramset}/contigs/{contig}/plots/metricsForPostImputationQC.sample.jpg"),
        temp("{paramset}/contigs/{contig}/plots/metricsForPostImputationQCChromosomeWide.sample.jpg"),
        temp("{paramset}/contigs/{contig}/plots/r2.goodonly.jpg"),
        temp(directory("{paramset}/contigs/{contig}/RData")),
        temp(directory("{paramset}/contigs/{contig}/input")),
        temp("{paramset}/contigs/{contig}/{contig}.vcf.gz"),
        tmpdir = temp(directory("{paramset}/contigs/{contig}/tmp"))
    log:
        stitch_log = "{paramset}/logs/{contig}.stitch.log",
        rename_log = "{paramset}/logs/{contig}.mv_stitchplots.log"
    params:
        chrom   = lambda wc: f"--chr={wc.contig}",
        start   = lambda wc: f"--regionStart={startpos}" if region else "",
        end     = lambda wc: f"--regionEnd={endpos}" if region else "",
        buffer  = lambda wc: f"--buffer={buffer}" if region else "",
        model   = lambda wc: f"--method={stitch_params[wc.paramset]['model']}",
        k       = lambda wc: f"--K={stitch_params[wc.paramset]['k']}",
        s       = lambda wc: f"--S={stitch_params[wc.paramset]['s']}",
        ngen    = lambda wc: f"--nGen={stitch_params[wc.paramset]['ngen']}",
        outdir  = lambda wc: "--outputdir=" + os.path.join(wc.paramset, "contigs", wc.contig),
        outfile = lambda wc: f"--output_filename={wc.contig}.vcf.gz",
        usebx   = lambda wc: f"--use_bx_tag={str(stitch_params[wc.paramset]['usebx']).upper()}",
        bxlimit = lambda wc: f"--bxTagUpperLimit={stitch_params[wc.paramset]['bxlimit']}",
        tmpdir  = lambda wc: "--tempdir=" + os.path.join(wc.paramset, "contigs", wc.contig, "tmp"),
        extra   = " ".join([f"{i}={j}" for i,j in extraparams.items()])
    threads:
        workflow.cores - 1
    conda:
        "envs/impute.yaml"
    container:
        f"docker://pdimens/harpy:impute_{VERSION}"        
    shell:
        """
        mkdir -p {output.tmpdir}
        STITCH.R --nCores={threads} --bamlist={input.bamlist} --posfile={input.infile} {params} 2> {log.stitch_log}
        {{
            cd {wildcards.paramset}/contigs/{wildcards.contig}/plots
            mv alphaMat.*all*.png alphaMat.all.png
            mv alphaMat.*normalized*.png alphaMat.normalized.png
            mv hapSum_log.*.png hapSum_log.png
            mv hapSum.*.png hapSum.png
            mv metricsForPostImputationQC.*sample.jpg metricsForPostImputationQC.sample.jpg
            mv metricsForPostImputationQCChromosomeWide*sample.jpg metricsForPostImputationQCChromosomeWide.sample.jpg
            mv r2*.goodonly.jpg r2.goodonly.jpg
        }} 2> {log.rename_log}
        """

rule index_vcf:
    input:
        "{paramset}/contigs/{contig}/{contig}.vcf.gz"
    output:
        "{paramset}/contigs/{contig}.vcf.gz.tbi",
        vcf   = "{paramset}/contigs/{contig}.vcf.gz",
        stats = "{paramset}/reports/data/contigs/{contig}.stats"
    shell:
        """
        cp {input} {output.vcf}
        tabix {output.vcf}
        bcftools stats -s "-" {input} > {output.stats}
        """

rule concat_list:
    input:
        cntg = collect("{{paramset}}/contigs/{contig}.vcf.gz", contig = contigs)
    output:
        bcf = temp("{paramset}/bcf.files")
    run:
        with open(output.bcf, "w") as fout:
            _ = fout.write("\n".join(input.cntg))

if len(contigs) == 1:
    rule bcf_conversion:
        input:
            collect("{{paramset}}/contigs/{contig}.vcf.gz.tbi", contig = contigs),
            vcf = collect("{{paramset}}/contigs/{contig}.vcf.gz", contig = contigs)
        output:
            "{paramset}/{paramset}.bcf.csi",
            bcf = "{paramset}/{paramset}.bcf"
        log:
            "logs/concat/{paramset}.concat.log"
        shell:
            "bcftools view -Ob --write-index -o {output.bcf} {input.vcf} 2> {log}"

else:
    rule merge_vcf:
        priority: 100
        input:
            collect("{{paramset}}/contigs/{contig}.vcf.gz.tbi", contig = contigs),
            files = "{paramset}/bcf.files"
        output:
            "{paramset}/{paramset}.bcf"
        log:
            "logs/concat/{paramset}.concat.log"
        threads:
            workflow.cores
        shell:
            "bcftools concat --threads {threads} -Ob -o {output} -f {input.files} 2> {log}"

    rule index_merged:
        input:
            "{paramset}/{paramset}.bcf"
        output:
            "{paramset}/{paramset}.bcf.csi"
        shell:
            "bcftools index {input}"

rule general_stats:
    input:
        "{paramset}/{paramset}.bcf.csi",
        bcf = "{paramset}/{paramset}.bcf"
    output:
        "{paramset}/reports/data/impute.stats"
    shell:
        "bcftools stats -s \"-\" {input.bcf} > {output}"

rule extract_region:
    input:
        "workflow/input/vcf/input.sorted.bcf.csi",
        orig = "workflow/input/vcf/input.sorted.bcf"
    output:
        temp("workflow/input/vcf/region.bcf.csi"),
        bcf = temp("workflow/input/vcf/region.bcf")
    params:
        f"-r {region}" if region else ""
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
    shell:
        """
        bcftools stats -s "-" {input.orig} {input.impute} | grep \"GCTs\" > {output.compare}
        bcftools query -f '%CHROM\\t%POS\\t%INFO/INFO_SCORE\\n' {input.impute} > {output.info_sc}
        """

rule contig_report:
    input:
        "{paramset}/contigs/{contig}/plots/alphaMat.all.png",
        "{paramset}/contigs/{contig}/plots/alphaMat.normalized.png",
        "{paramset}/contigs/{contig}/plots/hapSum_log.png",
        "{paramset}/contigs/{contig}/plots/hapSum.png",
        "{paramset}/contigs/{contig}/plots/metricsForPostImputationQC.sample.jpg",
        "{paramset}/contigs/{contig}/plots/metricsForPostImputationQCChromosomeWide.sample.jpg",
        "{paramset}/contigs/{contig}/plots/r2.goodonly.jpg",
        statsfile = "{paramset}/reports/data/contigs/{contig}.stats",
        ipynb = "workflow/stitch_collate.ipynb"
    output:
        tmp = temp("{paramset}/reports/{contig}.{paramset}.tmp.ipynb")
        ipynb = "{paramset}/reports/{contig}.{paramset}.ipynb"
    log:
        logfile = "{paramset}/logs/reports/{contig}.stitch.log"
    params:
        statsfile = lambda wc: "-p statsfile " + os.path.abspath("{wc.paramset}/reports/data/contigs/{wc.contig}.stats"),
        plotdir = lambda wc: "-p plotdir " + os.path.abspath(f"{wc.paramset}/contigs/{wc.contig}/plots"),
        model   = lambda wc: f"-p model {stitch_params[wc.paramset]['model']}",
        usebx   = lambda wc: f"-p usebx {stitch_params[wc.paramset]['usebx']}",
        bxlimit = lambda wc: f"-p bxlimit {stitch_params[wc.paramset]['bxlimit']}",
        k       = lambda wc: f"-p k {stitch_params[wc.paramset]['k']}",
        s       = lambda wc: f"-p s {stitch_params[wc.paramset]['s']}",
        ngen    = lambda wc: f"-p ngen {stitch_params[wc.paramset]['ngen']}",
        extra   = f"-P extra:{stitch_extra}"
    shell:
        """
        {{
            papermill -k python3 --no-progress-bar --log-level ERROR {input.ipynb} {output.tmp} {params}
            process-noteobok {wildcards.contig} {wildcards.paramset} {output.tmp}
        }} 2> {log} > {output.ipynb}
        """

rule impute_reports:
    input:
        comparison = "{paramset}/reports/data/impute.compare.stats",
        infoscore = "{paramset}/reports/data/impute.infoscore",
        ipynb = "workflow/impute.ipynb"
    output:
        tmp = temp("{paramset}/reports/{paramset}.summary.tmp.ipynb")
        ipynb = "{paramset}/reports/{paramset}.summary.ipynb"
    log:
        "{paramset}/logs/reports/imputestats.log"
    params:
        basedir = lambda wc: "-p basedir " + os.path.abspath("{wc.paramset}/reports/data/"),
        model   = lambda wc: f"-p model:{stitch_params[wc.paramset]['model']}",
        usebx   = lambda wc: f"-p usebx:{stitch_params[wc.paramset]['usebx']}",
        bxlimit = lambda wc: f"-p bxlimit:{stitch_params[wc.paramset]['bxlimit']}",
        k       = lambda wc: f"-p k:{stitch_params[wc.paramset]['k']}",
        s       = lambda wc: f"-p s:{stitch_params[wc.paramset]['s']}",
        ngen    = lambda wc: f"-p ngen:{stitch_params[wc.paramset]['ngen']}",
        extra   = f"-p extra:{stitch_extra}"
    shell:
        """
        {{
            papermill -k python3 --no-progress-bar --log-level ERROR {input.ipynb} {output.tmp} {params}
            process-noteobok {wildcards.paramset} {output.tmp}
        }} 2> {log} > {output.ipynb}
        """

rule all:
    default_target: True
    input: 
        vcf = collect("{paramset}/{paramset}.bcf", paramset = list(stitch_params.keys())),
        agg_report = collect("{paramset}/reports/{paramset}.summary.ipynb", paramset = stitch_params.keys()) if not skip_reports else [],
        contig_report = collect("{paramset}/reports/{contig}.{paramset}.ipynb", paramset = stitch_params.keys(), contig = contigs) if not skip_reports else []
