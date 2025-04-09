containerized: "docker://pdimens/harpy:latest"

import os
import re
import logging
from pathlib import Path

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake_log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+",
    population = r"[a-zA-Z0-9._-]+"

envdir       = os.path.join(os.getcwd(), "workflow", "envs")
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

rule preproc_groups:
    input:
        groupfile
    output:
        "workflow/sample.groups"
    run:
        with open(input[0], "r") as infile, open(output[0], "w") as outfile:
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
        temp("workflow/input/{population}.unsort.bam")
    log:
        "logs/concat_groups/{population}.concat.log"
    container:
        None
    shell:
        "concatenate_bam.py -b {input.bamlist} > {output} 2> {log}"

rule sort_groups:
    input:
        "workflow/input/{population}.unsort.bam"
    output:
        bam = ("workflow/input/{population}.bam"),
        bai = ("workflow/input/{population}.bam.bai")
    log:
        "logs/samtools/sort/{population}.sort.log"
    resources:
        mem_mb = 2000
    threads:
        10
    container:
        None
    shell:
        "samtools sort -@ {threads} -O bam -l 0 -m {resources.mem_mb}M --write-index -o {output.bam}##idx##{output.bai} {input} 2> {log}"

rule naibr_config:
    input:
        "workflow/input/{population}.bam"
    output:
        "workflow/config/{population}.naibr"
    params:
        lambda wc: wc.get("population"),
        min(10, workflow.cores - 1)
    run:
        with open(output[0], "w") as conf:
            _ = conf.write(f"bam_file={input[0]}\n")
            _ = conf.write(f"outdir={params[0]}\n")
            _ = conf.write(f"prefix={params[0]}\n")
            _ = conf.write(f"threads={params[1]}\n")
            for i in argdict:
                _ = conf.write(f"{i}={argdict[i]}\n")

rule call_variants:
    input:
        bam   = "workflow/input/{population}.bam",
        bai   = "workflow/input/{population}.bam.bai",
        conf  = "workflow/config/{population}.naibr"
    output:
        temp(directory("naibrlog")),
        bedpe = temp("{population}/{population}.bedpe"),
        refmt = temp("{population}/{population}.reformat.bedpe"),
        vcf   = temp("{population}/{population}.vcf"),
        log   = temp("{population}/{population}.log")
    log:
        "logs/naibr/{population}.naibr.log"
    threads:
        min(10, workflow.cores - 1)
    conda:
        f"{envdir}/variants.yaml"
    shell:
        "naibr {input.conf} > {log} 2>&1"

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
    container:
        None
    shell:
        """
        infer_sv.py {input.bedpe} -f {output.fail} > {output.bedpe}
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

rule process_genome:
    input:
        genomefile
    output: 
        workflow_geno
    container:
        None
    shell: 
        "seqtk seq {input} > {output}"

rule index_genome:
    input: 
        workflow_geno
    output: 
        f"{workflow_geno}.fai"
    log:
        f"{workflow_geno}.faidx.log"
    container:
        None
    shell:
        "samtools faidx --fai-idx {output} {input} 2> {log}"

rule report_config:
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
        f"{envdir}/r.yaml"
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        FAIDX=$(realpath {input.faidx})
        BEDPE=$(realpath {input.bedpe})
        quarto render {output.qmd} --log {log} --quiet -P faidx:$FAIDX -P bedpe:$BEDPE {params}
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
        f"{envdir}/r.yaml"
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        FAIDX=$(realpath {input.faidx})
        INPATH=$(realpath {params.bedpedir})
        quarto render {output.qmd} --log {log} --quiet -P faidx:$FAIDX -P bedpedir:$INPATH {params.contigs}
        """

rule workflow_summary:
    default_target: True
    input:
        bedpe = collect("bedpe/{pop}.bedpe", pop = populations),
        bedpe_agg = collect("{sv}.bedpe", sv = ["inversions", "deletions","duplications"]),
        reports = collect("reports/{pop}.naibr.html", pop = populations) if not skip_reports else [],
        agg_report = "reports/naibr.summary.html" if not skip_reports else []
    run:
        summary = ["The harpy sv naibr workflow ran using these parameters:"]
        summary.append(f"The provided reference genome: {bn}")
        concat = "The alignments were concatenated using:\n"
        concat += "\tconcatenate_bam.py -o groupname.bam -b samples.list"
        summary.append(concat)
        naibr = "naibr variant calling ran using these configurations:\n"
        naibr += "\tbam_file=BAMFILE\n"
        naibr += "\tprefix=PREFIX\n"
        naibr += "\toutdir=Variants/naibr/PREFIX\n"
        naibr += "\n\t".join([f"{k}={v}" for k,v in argdict.items()])
        summary.append(naibr)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['snakemake_command']}"
        summary.append(sm)
        with open("workflow/sv.naibr.summary", "w") as f:
            f.write("\n\n".join(summary))