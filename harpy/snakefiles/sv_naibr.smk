import os
import re
from pathlib import Path
from harpy.common.file_ops import naibr_extra, pop_manifest

localrules: all, naibr_config, aggregate_variants
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

WORKFLOW   = config.get('Workflow') or {}
PARAMETERS = config.get('Parameters') or {}
REPORTS    = WORKFLOW.get("reports") or {} 
INPUTS     = config['Inputs']
VERSION    = WORKFLOW.get('harpy-version', 'latest')

skip_reports = REPORTS.get("skip", False)
plot_contigs = REPORTS.get("plot-contigs", 'default')
extra        = PARAMETERS.get("extra", None) 
min_size     = PARAMETERS.get("min-size", 1000)
min_barcodes = PARAMETERS.get("min-barcodes", 2)
min_quality  = PARAMETERS.get("min-map-quality", 30)
mol_dist     = PARAMETERS.get("molecule-distance", 100000)
genomefile   = INPUTS["reference"]
bamlist      = INPUTS["alignments"]
# attempt to get processed, then source, then nothing
grp          = INPUTS.get("groupings") or {}
if grp:
    groupings = grp.get("processed", [])
    if not os.path.isfile(groupings):
        groupings.get("source") or []
else:
    groupings = []

plot_contigs  = ",".join(plot_contigs) if isinstance(plot_contigs, list) else plot_contigs
bamdict       = dict(zip(bamlist, bamlist))
popdict       = pop_manifest(groupings, bamlist) if groupings else None
populations   = popdict.keys() if groupings else None
target        = populations if groupings else {Path(i).stem for i in bamlist}
bn            = os.path.basename(genomefile)
workflow_geno = f"workflow/reference/{bn[:-3]}" if bn.lower().endswith(".gz") else f"workflow/reference/{bn}"
argdict = naibr_extra(
    {"min_mapq" : min_quality, "d" : mol_dist, "min_sv" : min_size, "k": min_barcodes},
    extra
)

def get_alignments(wildcards):
    """returns a list with the bam file for the sample based on wildcards.sample"""
    r = re.compile(fr".*/({wildcards.sample})\.(bam|sam)$", flags = re.IGNORECASE)
    aln = list(filter(r.match, bamlist))
    return aln[0]

def get_align_index(wildcards):
    """returns a list with the bai index file for the sample based on wildcards.sample"""
    r = re.compile(fr"(.*/{wildcards.sample})\.(bam|sam)$", flags = re.IGNORECASE)
    aln = list(filter(r.match, bamlist))
    return aln[0] + ".bai"

rule process_reference:
    input:
        genomefile
    output: 
        geno = workflow_geno,
        fai  = f"{workflow_geno}.fai"
    log:
        f"{workflow_geno}.preprocess.log"
    shell: 
        """
        {{
            XXseqtk seq {input} > {output.geno}
            samtools faidx --fai-idx {output.fai} {output.geno}
        }} 2> {log}
        """

if not groupings:
    rule index_alignments:
        input:
            lambda wc: bamdict[wc.bam]
        output:
            "{bam}.bai"
        shell:
            "samtools index {input}"

rule concat_groups:
    input: 
        bamfiles = lambda wc: collect("{samples}", samples = popdict[wc.sample]) 
    output:
        bam = temp("workflow/input/{sample}.bam"),
        bai = temp("workflow/input/{sample}.bam.bai")
    log:
        "logs/concat_groups/{sample}.concat.log"
    resources:
        mem_mb = 2000
    threads:
        workflow.cores
    shell:
        """
        {{
            djinn sam concat {input} | 
            samtools sort -@ {threads} -O bam -l 0 -m {resources.mem_mb}M --write-index -o {output.bam}##idx##{output.bai}
        }} 2> {log}
        """

rule naibr_config:
    input:
        bam = "workflow/input/{sample}.bam" if popdict else get_alignments
    output:
        cfg = "workflow/config/{sample}.naibr"
    params:
        popu = lambda wc: wc.get("sample"),
        thd = min(10, workflow.cores - 1)
    run:
        with open(output.cfg, "w") as conf:
            _ = conf.write(f"bam_file={input.bam}\n")
            _ = conf.write(f"outdir={params.popu}\n")
            _ = conf.write(f"prefix={params.popu}\n")
            _ = conf.write(f"threads={params.thd}\n")
            for i in argdict:
                _ = conf.write(f"{i}={argdict[i]}\n")

rule call_variants:
    input:
        bam   = "workflow/input/{sample}.bam" if popdict else get_alignments,
        bai   = "workflow/input/{sample}.bam.bai" if popdict else get_align_index,
        conf  = "workflow/config/{sample}.naibr"
    output:
        bedpe = temp("{sample}/{sample}.bedpe"),
        refmt = temp("{sample}/{sample}.reformat.bedpe"),
        vcf   = temp("{sample}/{sample}.vcf"),
        log   = temp("{sample}/{sample}.log")
    log:
        "logs/naibr/{sample}.naibr.log"
    threads:
        min(10, workflow.cores - 1)
    conda:
        "envs/variants.yaml"
    container:
        f"docker://pdimens/harpy:variants_{VERSION}"
    shell:
        "naibr {input.conf} > {log} 2>&1 && rm -rf naibrlog"

rule infer_variants:
    priority: 100
    input:
        bedpe = "{sample}/{sample}.bedpe",
        refmt = "{sample}/{sample}.reformat.bedpe",
        vcf   = "{sample}/{sample}.vcf"
    output:
        bedpe = "bedpe/{sample}.bedpe",
        refmt = "IGV/{sample}.reformat.bedpe",
        fail  = "bedpe/qc_fail/{sample}.fail.bedpe",
        vcf   = "vcf/{sample}.vcf" 
    shell:
        """
        harpy-utils infer-sv {input.bedpe} -f {output.fail} > {output.bedpe}
        cp {input.refmt} {output.refmt}
        cp {input.vcf} {output.vcf}
        """

rule aggregate_variants:
    input:
        collect("bedpe/{sample}.bedpe", sample = target)
    output:
        "inversions.bedpe",
        "deletions.bedpe",
        "duplications.bedpe"
    run:
        with open(output[0], "w") as inversions, open(output[1], "w") as deletions, open(output[2], "w") as duplications:
            header = ["Population","Chr1","Break1","Chr2","Break2","SplitMolecules","DiscordantReads","Orientation","Haplotype","Score","PassFilter","SV"]
            _ = inversions.write("\t".join(header) + "\n")
            _ = deletions.write("\t".join(header) + "\n")
            _ = duplications.write("\t".join(header) + "\n")
            for varfile in input:
                samplename = Path(varfile).stem
                with open(varfile, "r") as f_in:
                    # read the header to skip it
                    f_in.readline()
                    # read the rest of it
                    while True:
                        line = f_in.readline()
                        if not line:
                            break
                        record = line.rstrip().split("\t")
                        if record[-1] == "inversion":
                            _ = inversions.write(f"{samplename}\t{line}")
                        elif record[-1] == "deletion":
                            _ = deletions.write(f"{samplename}\t{line}")
                        elif record[-1] == "duplication":
                            _ = duplications.write(f"{samplename}\t{line}")

rule report:
    input: 
        faidx = f"{workflow_geno}.fai",
        stats = collect("{var}.bedpe", var = ['inversions', 'deletions', 'duplications']),
        ipynb = "workflow/sv.ipynb"
    output:
        tmp = temp("reports/naibr.summary.tmp.ipynb"),
        ipynb = "reports/naibr.summary.ipynb"
    log:
        "logs/report.log"
    params:
        f"-p indir {os.getcwd()}",
        f"-p faidx " + os.path.abspath(f"{workflow_geno}.fai"),
        f"-p contigs {plot_contigs}" if plot_contigs != "default" else ""
    shell:
        """
        {{
            papermill -k xpython --no-progress-bar --log-level ERROR {input.ipynb} {output.tmp} {params}
            harpy-utils process-notebook {output.tmp} NAIBR
        }} 2> {log} > {output.ipynb}
        """

rule all:
    default_target: True
    input:
        bedpe = collect("bedpe/{sample}.bedpe", sample = target),
        bedpe_agg = collect("{sv}.bedpe", sv = ["inversions", "deletions","duplications"]),
        agg_report = "reports/naibr.summary.ipynb" if not skip_reports else []
