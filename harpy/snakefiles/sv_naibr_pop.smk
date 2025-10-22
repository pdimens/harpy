import os
import re
import logging
from pathlib import Path

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake"]["log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+",
    population = r"[a-zA-Z0-9._-]+"

genomefile   = config["inputs"]["reference"]
bamlist      = config["inputs"]["alignments"]
groupfile    = config["inputs"]["groupings"]
extra        = config.get("extra", None) 
min_size       = config["min_size"]
min_barcodes = config["min_barcodes"]
min_quality  = config["min_quality"]
mol_dist     = config["molecule_distance"]
skip_reports  = config["reports"]["skip"]
plot_contigs = config["reports"]["plot_contigs"]    
bn           = os.path.basename(genomefile)
if bn.lower().endswith(".gz"):
    workflow_geno = f"workflow/reference/{bn[:-3]}"
else:
    workflow_geno = f"workflow/reference/{bn}"

def process_args(args):
    argsDict = {
        "min_mapq" : min_quality,
        "d"        : mol_dist,
        "min_sv"   : min_size,
        "k"        : min_barcodes
    }
    if args:
        words = [i for i in re.split(r"\s|=", args) if len(i) > 0]
        for i in zip(words[::2], words[1::2]):
            if "blacklist" in i or "candidates" in i:
                argsDict[i[0].lstrip("-")] = i[1]
    return argsDict

argdict = process_args(extra)

# create dictionary of population => filenames
## this makes it easier to set the snakemake rules/wildcards
def pop_manifest(groupingfile, filelist):
    d = {}
    with open(groupingfile) as f:
        for line in f:
            samp, pop = line.rstrip().split()
            if samp.lstrip().startswith("#"):
                continue
            r = re.compile(fr".*/({samp.lstrip()})\.(bam|sam)$", flags = re.IGNORECASE)
            sampl = list(filter(r.match, filelist))[0]
            if pop not in d.keys():
                d[pop] = [sampl]
            else:
                d[pop].append(sampl)
    return d

popdict = pop_manifest(groupfile, bamlist)
populations = popdict.keys()

rule preprocess_groups:
    input:
        grp = groupfile
    output:
        grp = "workflow/sample.groups"
    run:
        with open(input.grp, "r") as infile, open(output.grp, "w") as outfile:
            _ = [outfile.write(i) for i in infile.readlines() if not i.lstrip().startswith("#")]

rule concat_list:
    input:
        "workflow/sample.groups"
    output:
        "workflow/merge_samples/{population}.list"
    run:
        with open(output[0], "w") as fout:
            for bamfile in popdict[wildcards.population]:
                _ = fout.write(bamfile + "\n")

rule concat_groups:
    input: 
        bamlist  = "workflow/merge_samples/{population}.list",
        bamfiles = lambda wc: collect("{sample}", sample = popdict[wc.population]) 
    output:
        bam = "workflow/input/{population}.bam",
        bai = "workflow/input/{population}.bam.bai"
    log:
        "logs/concat_groups/{population}.concat.log"
    resources:
        mem_mb = 2000
    threads:
        10
    shell:
        """
        {{
            concatenate_bam -b {input.bamlist} |
            samtools sort -@ {threads} -O bam -l 0 -m {resources.mem_mb}M --write-index -o {output.bam}##idx##{output.bai}
        }} 2> {log}
        """

rule naibr_config:
    input:
        bam = "workflow/input/{population}.bam"
    output:
        cfg = "workflow/config/{population}.naibr"
    params:
        popu = lambda wc: wc.get("population"),
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
        bam   = "workflow/input/{population}.bam",
        bai   = "workflow/input/{population}.bam.bai",
        conf  = "workflow/config/{population}.naibr"
    output:
        bedpe = temp("{population}/{population}.bedpe"),
        refmt = temp("{population}/{population}.reformat.bedpe"),
        vcf   = temp("{population}/{population}.vcf"),
        log   = temp("{population}/{population}.log")
    log:
        "logs/naibr/{population}.naibr.log"
    threads:
        min(10, workflow.cores - 1)
    conda:
        "envs/variants.yaml"
    container:
        "docker://pdimens/harpy:variants_latest"
    shell:
        "naibr {input.conf} > {log} 2>&1 && rm -rf naibrlog"

rule infer_variants:
    priority: 100
    input:
        bedpe = "{population}/{population}.bedpe",
        refmt = "{population}/{population}.reformat.bedpe",
        vcf   = "{population}/{population}.vcf"
    output:
        bedpe = "bedpe/{population}.bedpe",
        refmt = "IGV/{population}.reformat.bedpe",
        fail  = "bedpe/qc_fail/{population}.fail.bedpe",
        vcf   = "vcf/{population}.vcf" 
    shell:
        """
        infer_sv {input.bedpe} -f {output.fail} > {output.bedpe}
        cp {input.refmt} {output.refmt}
        cp {input.vcf} {output.vcf}
        """

rule aggregate_variants:
    input:
        collect("bedpe/{population}.bedpe", population = populations)
    output:
        "inversions.bedpe",
        "deletions.bedpe",
        "duplications.bedpe"
    run:
        from pathlib import Path
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

rule preprocess_reference:
    input:
        genomefile
    output: 
        geno = workflow_geno,
        fai = f"{workflow_geno}.fai"
    log:
        f"{workflow_geno}.preprocess.log"
    shell: 
        """
        {{
            seqtk seq {input} > {output.geno}
            samtools faidx --fai-idx {output.fai} {output.geno}
        }} 2> {log}
        """

rule configure_report:
    input:
        yaml = "workflow/report/_quarto.yml",
        scss = "workflow/report/_harpy.scss"
    output:
        yaml = temp("reports/_quarto.yml"),
        scss = temp("reports/_harpy.scss")
    run:
        import shutil
        for i,o in zip(input,output):
            shutil.copy(i,o)

rule group_reports:
    input: 
        "reports/_quarto.yml",
        "reports/_harpy.scss",
        faidx = f"{workflow_geno}.fai",
        bedpe = "bedpe/{population}.bedpe",
        qmd   = "workflow/report/naibr.qmd"
    output:
        report = "reports/{population}.naibr.html",
        qmd = temp("reports/{population}.naibr.qmd")
    log:
        "logs/reports/{population}.report.log"
    params:
        sample= lambda wc: "-P sample:" + wc.get('population'),
        contigs= f"-P contigs:{plot_contigs}"
    conda:
        "envs/report.yaml"
    container:
        "docker://pdimens/harpy:report_latest"
    retries:
        3
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        FAIDX=$(realpath {input.faidx})
        BEDPE=$(realpath {input.bedpe})
        quarto render {output.qmd} --no-cache --log {log} --quiet -P faidx:$FAIDX -P bedpe:$BEDPE {params}
        """

rule aggregate_report:
    input: 
        "reports/_quarto.yml",
        "reports/_harpy.scss",
        faidx = f"{workflow_geno}.fai",
        bedpe = collect("bedpe/{pop}.bedpe", pop = populations),
        qmd   = "workflow/report/naibr_pop.qmd"
    output:
        report = "reports/naibr.summary.html",
        qmd = temp("reports/naibr.summary.qmd")
    log:
        "logs/reports/summary.report.log"
    params:
        bedpedir = "bedpe",
        contigs = f"-P contigs:{plot_contigs}"
    conda:
        "envs/report.yaml"
    container:
        "docker://pdimens/harpy:report_latest"
    retries:
        3
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        FAIDX=$(realpath {input.faidx})
        INPATH=$(realpath {params.bedpedir})
        quarto render {output.qmd} --no-cache --log {log} --quiet -P faidx:$FAIDX -P bedpedir:$INPATH {params.contigs}
        """

rule all:
    default_target: True
    input:
        bedpe = collect("bedpe/{pop}.bedpe", pop = populations),
        bedpe_agg = collect("{sv}.bedpe", sv = ["inversions", "deletions","duplications"]),
        reports = collect("reports/{pop}.naibr.html", pop = populations) if not skip_reports else [],
        agg_report = "reports/naibr.summary.html" if not skip_reports else []
