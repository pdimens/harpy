containerized: "docker://pdimens/harpy:latest"

import os
import re
import sys
import multiprocessing
import logging as pylogging
from pathlib import Path

envdir       = os.getcwd() + "/.harpy_envs"
genomefile   = config["inputs"]["genome"]
bamlist      = config["inputs"]["alignments"]
groupfile    = config["inputs"]["groupings"]
extra        = config.get("extra", None) 
min_sv       = config["min_sv"]
min_barcodes = config["min_barcodes"]
min_quality  = config["min_quality"]
mol_dist     = config["molecule_distance"]
outdir       = config["output_directory"]
skipreports  = config["skip_reports"]
snakemake_log = config["snakemake_log"]
bn           = os.path.basename(genomefile)
if bn.lower().endswith(".gz"):
    bn = bn[:-3]

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+",
    population = "[a-zA-Z0-9._-]+"

onstart:
    extra_logfile_handler = pylogging.FileHandler(snakemake_log)
    logger.logger.addHandler(extra_logfile_handler)

def process_args(args):
    argsDict = {
        "min_mapq" : min_quality,
        "d"        : mol_dist,
        "min_sv"   : min_sv,
        "k"        : min_barcodes
    }
    if args:
        words = [i for i in re.split(r"\s|=", args) if len(i) > 0]
        for i in zip(words[::2], words[1::2]):
            if "blacklist" in i or "candidates" in i:
                argsDict[i[0].lstrip("-")] = i[1]
    return argsDict

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
        outdir + "/workflow/sample.groups"
    run:
        with open(input[0], "r") as infile, open(output[0], "w") as outfile:
            _ = [outfile.write(i) for i in infile.readlines() if not i.lstrip().startswith("#")]

rule concat_list:
    input:
        outdir + "/workflow/sample.groups"
    output:
        outdir + "/workflow/merge_samples/{population}.list"
    run:
        with open(output[0], "w") as fout:
            for bamfile in popdict[wildcards.population]:
                _ = fout.write(bamfile + "\n")

rule concat_groups:
    input: 
        bamlist  = outdir + "/workflow/merge_samples/{population}.list",
        bamfiles = lambda wc: collect("{sample}", sample = popdict[wc.population]) 
    output:
        temp(outdir + "/workflow/input/{population}.unsort.bam")
    log:
        outdir + "/logs/concat_groups/{population}.concat.log"
    threads:
        1
    container:
        None
    shell:
        "concatenate_bam.py -o {output} -b {input.bamlist} 2> {log}"

rule sort_groups:
    input:
        outdir + "/workflow/input/{population}.unsort.bam"
    output:
        bam = (outdir + "/workflow/input/{population}.bam"),
        bai = (outdir + "/workflow/input/{population}.bam.bai")
    log:
        outdir + "/logs/samtools/sort/{population}.sort.log"
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
        outdir + "/workflow/input/{population}.bam"
    output:
        outdir + "/workflow/config/{population}.naibr"
    params:
        lambda wc: wc.get("population"),
        min(10, workflow.cores - 1)
    run:
        argdict = process_args(extra)
        with open(output[0], "w") as conf:
            _ = conf.write(f"bam_file={input[0]}\n")
            _ = conf.write(f"outdir={outdir}/{params[0]}\n")
            _ = conf.write(f"prefix={params[0]}\n")
            _ = conf.write(f"threads={params[1]}\n")
            for i in argdict:
                _ = conf.write(f"{i}={argdict[i]}\n")

rule call_variants:
    input:
        bam   = outdir + "/workflow/input/{population}.bam",
        bai   = outdir + "/workflow/input/{population}.bam.bai",
        conf  = outdir + "/workflow/config/{population}.naibr"
    output:
        bedpe = outdir + "/{population}/{population}.bedpe",
        refmt = outdir + "/{population}/{population}.reformat.bedpe",
        vcf   = outdir + "/{population}/{population}.vcf"
    log:
        outdir + "/logs/naibr/{population}.naibr.log"
    threads:
        min(10, workflow.cores - 1)
    conda:
        f"{envdir}/sv.yaml"
    shell:
        "naibr {input.conf} > {log} 2>&1"

rule infer_variants:
    input:
        bedpe = outdir + "/{population}/{population}.bedpe",
        refmt = outdir + "/{population}/{population}.reformat.bedpe",
        vcf   = outdir + "/{population}/{population}.vcf"
    output:
        bedpe = outdir + "/bedpe/{population}.bedpe",
        refmt = outdir + "/IGV/{population}.reformat.bedpe",
        fail  = outdir + "/bedpe/qc_fail/{population}.fail.bedpe",
        vcf   = outdir + "/vcf/{population}.vcf" 
    params:
        outdir = lambda wc: outdir + "/" + wc.get("population")
    container:
        None
    shell:
        """
        infer_sv.py {input.bedpe} -f {output.fail} > {output.bedpe}
        mv {input.refmt} {output.refmt} &&
        mv {input.vcf} {output.vcf} &&
        rm -rf {params.outdir}
        """

rule aggregate_variants_variants:
    input:
        collect(outdir + "/bedpe/{population}.bedpe", population = populations)
    output:
        outdir + "/inversions.bedpe",
        outdir + "/deletions.bedpe",
        outdir + "/duplications.bedpe"
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

rule setup_genome:
    input:
        genomefile
    output: 
        f"Genome/{bn}"
    container:
        None
    shell: 
        "seqtk seq {input} > {output}"

rule faidx_genome:
    input: 
        f"Genome/{bn}"
    output: 
        f"Genome/{bn}.fai"
    log:
        f"Genome/{bn}.faidx.log"
    container:
        None
    shell:
        "samtools faidx --fai-idx {output} {input} 2> {log}"

rule group_reports:
    input:
        fai   = f"Genome/{bn}.fai",
        bedpe = outdir + "/bedpe/{population}.bedpe"
    output:
        outdir + "/reports/{population}.naibr.html"
    log:
        logfile = outdir + "/logs/reports/{population}.report.log"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/naibr.Rmd"

rule aggregate_report:
    input:
        fai   = f"Genome/{bn}.fai",
        bedpe = collect(outdir + "/bedpe/{pop}.bedpe", pop = populations)
    output:
        outdir + "/reports/naibr.pop.summary.html"
    log:
        logfile = outdir + "/logs/reports/summary.report.log"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/naibr_pop.Rmd"

rule workflow_summary:
    default_target: True
    input:
        bedpe = collect(outdir + "/bedpe/{pop}.bedpe", pop = populations),
        bedpe_agg = collect(outdir + "/{sv}.bedpe", sv = ["inversions", "deletions","duplications"]),
        reports = collect(outdir + "/reports/{pop}.naibr.html", pop = populations) if not skipreports else [],
        agg_report = outdir + "/reports/naibr.pop.summary.html" if not skipreports else []
    run:
        os.system(f"rm -rf {outdir}/naibrlog")
        argdict = process_args(extra)
        with open(outdir + "/workflow/sv.naibr.summary", "w") as f:
            _ = f.write("The harpy sv naibr workflow ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write(f"The sample grouping file: {groupfile}\n\n")
            _ = f.write("The alignments were concatenated using:\n")
            _ = f.write("    concatenate_bam.py -o groupname.bam -b samples.list\n")
            _ = f.write("naibr variant calling ran using these configurations:\n")
            _ = f.write(f"    bam_file=BAMFILE\n")
            _ = f.write(f"    prefix=PREFIX\n")
            _ = f.write(f"    outdir=Variants/naibr/PREFIX\n")
            for i in argdict:
                _ = f.write(f"    {i}={argdict[i]}\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")
