import os
import re

localrules: all, alignment_list
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+",
    paramset = r"[^/]+",
    contig = r"[^/]+"

WORKFLOW   = config.get('Workflow') or {}
PARAMETERS = config.get('Parameters') or {}
REPORTS    = WORKFLOW.get("reports") or {} 
INPUTS     = config['Inputs']
VERSION    = WORKFLOW.get('harpy-version', 'latest')

skip_reports  = REPORTS.get("skip", False)
region        = PARAMETERS.get("region", None)
window        = PARAMETERS.get("window-size", None)
buffer        = PARAMETERS.get("buffer", 100000)
stitch_params = PARAMETERS["stitch"]
stitch_extra  = PARAMETERS.get("extra", "None")
grid_size     = PARAMETERS.get("grid-size", 1)
variantfile   = INPUTS["vcf"]
bamlist       = INPUTS["alignments"]

def makewindows(_c_len, windowsize):
    """create vectors of specified windows"""
    end = min(_c_len, windowsize)
    starts = [1]
    ends = [end]
    while end < _c_len:
        end = min(end + windowsize, _c_len)
        ends.append(end)
        starts.append(starts[-1] + windowsize)
    return [f"{i}-{j}" for i,j in zip(starts,ends)]

class Contig:
    def __init__(self, name, end, start = None, window = None):
        self.name: str = name
        if start:
            self.regions = [f"{start}-{end}"]
        elif window:
            self.regions = makewindows(end, window)
        else:
            self.regions = [f"1-{end}"]

def extract_contig_info():
    for contig in contigs.values():
        for region in contig.regions:
            yield {"contig": contig.name, "region": region}

bamdict = dict(zip(bamlist, bamlist))
contigs: dict[str, Contig] = {}
if region:
    cntg = region.split(":")[0]
    start,end = region.split('-')
    contigs[cntg] = Contig(cntg, int(end), start = int(start))
else:
    with open(INPUTS["biallelic-contigs"], "r") as f:
        for line in f:
            cntg,bp = line.strip().split()
            contigs[cntg] = Contig(cntg, int(bp), window = window)

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
        temp("{paramset}/contigs/{contig}/{region}/plots/alphaMat.all.png"),
        temp("{paramset}/contigs/{contig}/{region}/plots/alphaMat.normalized.png"),
        temp("{paramset}/contigs/{contig}/{region}/plots/hapSum_log.png"),
        temp("{paramset}/contigs/{contig}/{region}/plots/hapSum.png"),
        temp("{paramset}/contigs/{contig}/{region}/plots/metricsForPostImputationQC.sample.jpg"),
        temp("{paramset}/contigs/{contig}/{region}/plots/metricsForPostImputationQCChromosomeWide.sample.jpg"),
        temp("{paramset}/contigs/{contig}/{region}/plots/r2.goodonly.jpg"),
        temp(directory("{paramset}/contigs/{contig}/{region}/RData")),
        temp(directory("{paramset}/contigs/{contig}/{region}/input")),
        temp("{paramset}/contigs/{contig}/{region}/{contig}.{region}.vcf.gz.tbi")
        vcf = temp("{paramset}/contigs/{contig}/{region}/{contig}.{region}.vcf.gz"),
        tmpdir = temp(directory("{paramset}/contigs/{contig}/{region}/tmp"))
    log:
        "{paramset}/logs/{contig}.{region}.stitch.log",
    params:
        model   = lambda wc: f"--method={stitch_params[wc.paramset]['model']}",
        k       = lambda wc: f"--K={stitch_params[wc.paramset]['k']}",
        s       = lambda wc: f"--S={stitch_params[wc.paramset]['s']}",
        ngen    = lambda wc: f"--nGen={stitch_params[wc.paramset]['ngen']}",
        outdir  = lambda wc: "--outputdir=" + os.path.join(wc.paramset, "contigs", wc.contig, wc.region),
        outfile = lambda wc: f"--output_filename={wc.contig}.{wc.region}.vcf.gz",
        usebx   = lambda wc: f"--use_bx_tag={str(stitch_params[wc.paramset]['usebx']).upper()}",
        bxlimit = lambda wc: f"--bxTagUpperLimit={stitch_params[wc.paramset]['bxlimit']}",
        tmpdir  = lambda wc: f"--tempdir=" + os.path.join(wc.paramset, "contigs", wc.contig, wc.region, "tmp"),
        extra   = " ".join([f"{i}={j}" for i,j in extraparams.items()]),
        chrom   = lambda wc: f"--chr={wc.contig}",
        start   = lambda wc: f"--regionStart={wc.region.split('-')[0]}" if region or window else "",
        end     = lambda wc: f"--regionEnd={wc.region.split('-')[1]}" if region or window else "",
        buffer  = lambda wc: f"--buffer={buffer}" if region or window else "",
    threads:
        workflow.cores - 1
    conda:
        "envs/impute.yaml"
    container:
        f"docker://pdimens/harpy:impute_{VERSION}"        
    shell:
        """
        {{
            mkdir -p {output.tmpdir}
            STITCH.R --nCores={threads} --bamlist={input.bamlist} --posfile={input.infile} {params}
            cd {wildcards.paramset}/contigs/{wildcards.contig}/{wildcards.region}/plots
            mv alphaMat.*all*.png alphaMat.all.png
            mv alphaMat.*normalized*.png alphaMat.normalized.png
            mv hapSum_log.*.png hapSum_log.png
            mv hapSum.*.png hapSum.png
            mv metricsForPostImputationQC.*sample.jpg metricsForPostImputationQC.sample.jpg
            mv metricsForPostImputationQCChromosomeWide*sample.jpg metricsForPostImputationQCChromosomeWide.sample.jpg
            mv r2*.goodonly.jpg r2.goodonly.jpg
            tabix {output.vcf}
        }} 2> {log}
        """

rule merge_regions:
    priority: 100
    input:
        collect("{{paramset}}/contigs/{contig}/{region}/{contig}.{region}.vcf.gz", zip,
            contig = [name for name, c in contigs.items() for _ in c.regions],
            region=[r for c in contigs.values() for r in c.regions]),
        collect("{{paramset}}/contigs/{contig}/{region}/{contig}.{region}.vcf.gz.tbi", zip,
            contig = [name for name, c in contigs.items() for _ in c.regions],
            region=[r for c in contigs.values() for r in c.regions]),
    output:
        "{paramset}/{paramset}.bcf.csi",
        bcf = "{paramset}/{paramset}.bcf",
        filelist = temp("{paramset}/bcf.files")
    log:
        "logs/concat/{paramset}.concat.log"
    threads:
        workflow.cores
    shell:
        """
        {{
            find {wildcards.paramset}/contigs -name '*.vcf.gz' > {output.filelist}
            bcftools concat -a --threads {threads} -f {output.filelist} |
            bcftools sort -Ob -o {output.bcf} --write-index
        }}  2> {log}
        """
        # --allow-overlaps

rule extract_original_region:
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

#rule contig_report:
#    input:
#        "{paramset}/contigs/{contig}.{region}/plots/alphaMat.all.png",
#        "{paramset}/contigs/{contig}.{region}/plots/alphaMat.normalized.png",
#        "{paramset}/contigs/{contig}.{region}/plots/hapSum_log.png",
#        "{paramset}/contigs/{contig}.{region}/plots/hapSum.png",
#        "{paramset}/contigs/{contig}.{region}/plots/metricsForPostImputationQC.sample.jpg",
#        "{paramset}/contigs/{contig}.{region}/plots/metricsForPostImputationQCChromosomeWide.sample.jpg",
#        "{paramset}/contigs/{contig}.{region}/plots/r2.goodonly.jpg",
#        "{paramset}/contigs/{contig}.{region}.vcf.gz.tbi",
#        vcf   = "{paramset}/contigs/{contig}.{region}.vcf.gz",
#        ipynb = "workflow/stitch_collate.ipynb"
#    output:
#        stats = "{paramset}/reports/data/contigs/{contig}.stats",
#        tmp = temp("{paramset}/reports/{contig}.{paramset}.tmp.ipynb"),
#        ipynb = "{paramset}/reports/{contig}.{paramset}.ipynb"
#    log:
#        logfile = "{paramset}/logs/reports/{contig}.stitch.log"
#    params:
#        stats   = lambda wc: "-p statsfile " + os.path.abspath(f"{wc.paramset}/reports/data/contigs/{wc.contig}.stats"),
#        plotdir = lambda wc: "-p plotdir " + os.path.abspath(f"{wc.paramset}/contigs/{wc.contig}/plots"),
#        model   = lambda wc: f"-p model {stitch_params[wc.paramset]['model']}",
#        usebx   = lambda wc: f"-p usebx {stitch_params[wc.paramset]['usebx']}",
#        bxlimit = lambda wc: f"-p bxlimit {stitch_params[wc.paramset]['bxlimit']}",
#        k       = lambda wc: f"-p k {stitch_params[wc.paramset]['k']}",
#        s       = lambda wc: f"-p s {stitch_params[wc.paramset]['s']}",
#        ngen    = lambda wc: f"-p ngen {stitch_params[wc.paramset]['ngen']}",
#        extra   = f"-p extra {stitch_extra}"
#    shell:
#        """
#        {{
#            bcftools stats -s "-" {input.vcf} > {output.stats}
#            papermill -k xpython --no-progress-bar --log-level ERROR {input.ipynb} {output.tmp} {params}
#            harpy-utils process-notebook {output.tmp} {wildcards.contig} {wildcards.paramset}
#        }} 2> {log} > {output.ipynb}
#        """

rule impute_reports:
    input:
        orig    = "workflow/input/vcf/input.sorted.bcf" if not region else "workflow/input/vcf/region.bcf",
        origidx = "workflow/input/vcf/input.sorted.bcf.csi" if not region else "workflow/input/vcf/region.bcf.csi",
        impute  = "{paramset}/{paramset}.bcf",
        idx     = "{paramset}/{paramset}.bcf.csi",
        ipynb = "workflow/impute.ipynb"
    output:
        comparison = "{paramset}/reports/data/impute.compare.stats",
        infoscore = temp("{paramset}/reports/data/impute.infoscore"),
        tmp = temp("{paramset}/reports/{paramset}.summary.tmp.ipynb"),
        ipynb = "{paramset}/reports/{paramset}.summary.ipynb"
    log:
        "{paramset}/logs/reports/imputestats.log"
    params:
        basedir = lambda wc: "-p basedir " + os.path.abspath(f"{wc.paramset}/reports/data/"),
        model   = lambda wc: f"-p model {stitch_params[wc.paramset]['model']}",
        usebx   = lambda wc: f"-p usebx {stitch_params[wc.paramset]['usebx']}",
        bxlimit = lambda wc: f"-p bxlimit {stitch_params[wc.paramset]['bxlimit']}",
        k       = lambda wc: f"-p k {stitch_params[wc.paramset]['k']}",
        s       = lambda wc: f"-p s {stitch_params[wc.paramset]['s']}",
        ngen    = lambda wc: f"-p ngen {stitch_params[wc.paramset]['ngen']}",
        extra   = f"-p extra {stitch_extra}"
    shell:
        """
        {{
            bcftools stats -s "-" {input.orig} {input.impute} | grep \"GCTs\" > {output.comparison}
            bcftools query -f '%CHROM\\t%POS\\t%INFO/INFO_SCORE\\n' {input.impute} > {output.infoscore}
            papermill -k xpython --no-progress-bar --log-level ERROR {input.ipynb} {output.tmp} {params}
            harpy-utils process-notebook {output.tmp} {wildcards.paramset}
        }} 2> {log} > {output.ipynb}
        """

rule all:
    default_target: True
    input: 
        vcf = collect("{paramset}/{paramset}.bcf", paramset = list(stitch_params.keys())),
        agg_report = collect("{paramset}/reports/{paramset}.summary.ipynb", paramset = stitch_params.keys()) if not skip_reports else [],
        #contig_report = collect("{paramset}/reports/{contig}.{paramset}.ipynb", paramset = stitch_params.keys(), contig = contigs) if not skip_reports else []
