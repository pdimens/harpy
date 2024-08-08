containerized: "docker://pdimens/harpy:latest"

import os
import re
import sys
import glob
import multiprocessing
import logging as pylogging
from datetime import datetime
from pathlib import Path
from rich import print as rprint
from rich.panel import Panel

envdir      = os.getcwd() + "/.harpy_envs"
genomefile  = config["inputs"]["genome"]
bamlist     = config["inputs"]["alignments"]
samplenames = {Path(i).stem for i in bamlist}
extra       = config.get("extra", None) 
mol_dist    = config["molecule_distance"]
min_sv      = config["min_sv"]
min_barcodes = config["min_barcodes"]
min_quality  = config["min_quality"]
outdir      = config["output_directory"]
skipreports = config["skip_reports"]
bn          = os.path.basename(genomefile)
genome_zip  = True if bn.lower().endswith(".gz") else False
bn_idx      = f"{bn}.gzi" if genome_zip else f"{bn}.fai"

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

## the log file ##
attempts = glob.glob(f"{outdir}/logs/snakemake/*.snakelog")
if not attempts:
    logfile = f"{outdir}/logs/snakemake/sv_naibr.run1." + datetime.now().strftime("%d_%m_%Y") + ".snakelog"
else:
    increment = sorted([int(i.split(".")[1].replace("run","")) for i in attempts])[-1] + 1
    logfile = f"{outdir}/logs/snakemake/sv_naibr.run{increment}." + datetime.now().strftime("%d_%m_%Y") + ".snakelog"

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

onstart:
    os.makedirs(f"{outdir}/logs/snakemake", exist_ok = True)
    extra_logfile_handler = pylogging.FileHandler(logfile)
    logger.logger.addHandler(extra_logfile_handler)

onerror:
    print("")
    rprint(
        Panel(
            f"The workflow has terminated due to an error. See the log file for more details:\n[bold]{outdir}/logs/snakemake/{dt_string}.snakelog[/bold]",
            title = "[bold]harpy sv naibr",
            title_align = "left",
            border_style = "red"
            ),
        file = sys.stderr
    )

onsuccess:
    print("")
    rprint(
        Panel(
            f"The workflow has finished successfully! Find the results in [bold]{outdir}/[/bold]",
            title = "[bold]harpy sv naibr",
            title_align = "left",
            border_style = "green"
            ),
        file = sys.stderr
    )

def sam_index(infile):
    """Use Samtools to index an input file, adding .bai to the end of the name"""
    if not os.path.exists(f"{infile}.bai"):
        subprocess.run(f"samtools index {infile} {infile}.bai".split())

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

rule index_alignments:
    input:
        bamlist
    output:
        [f"{i}.bai" for i in bamlist]
    threads:
        workflow.cores
    message:
        "Indexing alignment files"
    run:
        with multiprocessing.Pool(processes=threads) as pool:
            pool.map(sam_index, input)

rule create_config:
    input:
        get_alignments
    output:
        outdir + "/workflow/input/{sample}.naibr"
    params:
        lambda wc: wc.get("sample"),
        min(10, workflow.cores)
    message:
        "Creating naibr config file: {wildcards.sample}"
    run:
        argdict = process_args(extra)
        with open(output[0], "w") as conf:
            _ = conf.write(f"bam_file={input[0]}\n")
            _ = conf.write(f"prefix={params[0]}\n")
            _ = conf.write(f"outdir={outdir}/{params[0]}\n")
            _ = conf.write(f"threads={params[1]}\n")
            for i in argdict:
                _ = conf.write(f"{i}={argdict[i]}\n")

rule call_sv:
    input:
        bam   = get_alignments,
        bai   = get_align_index,
        conf  = outdir + "/workflow/input/{sample}.naibr"
    output:
        bedpe = outdir + "/{sample}/{sample}.bedpe",
        refmt = outdir + "/{sample}/{sample}.reformat.bedpe",
        vcf   = outdir + "/{sample}/{sample}.vcf"
    log:
        outdir + "/logs/{sample}.naibr.log"
    threads:
        10
    conda:
        f"{envdir}/sv.yaml"     
    message:
        "Calling variants: {wildcards.sample}"
    shell:
        "naibr {input.conf} > {log} 2>&1"

rule infer_sv:
    input:
        bedpe = outdir + "/{sample}/{sample}.bedpe",
        refmt = outdir + "/{sample}/{sample}.reformat.bedpe",
        vcf   = outdir + "/{sample}/{sample}.vcf"
    output:
        bedpe = outdir + "/bedpe/{sample}.bedpe",
        refmt = outdir + "/IGV/{sample}.reformat.bedpe",
        fail  = outdir + "/bedpe/qc_fail/{sample}.fail.bedpe",
        vcf   = outdir + "/vcf/{sample}.vcf" 
    params:
        outdir = lambda wc: outdir + "/" + wc.get("sample")
    container:
        None
    message:
        "Inferring variants from naibr output: {wildcards.sample}"
    shell:
        """
        infer_sv.py {input.bedpe} -f {output.fail} > {output.bedpe}
        mv {input.refmt} {output.refmt} &&
        mv {input.vcf} {output.vcf} &&
        rm -rf {params.outdir}
        """

rule merge_variants:
    input:
        collect(outdir + "/bedpe/{sample}.bedpe", sample = samplenames)
    output:
        outdir + "/inversions.bedpe",
        outdir + "/deletions.bedpe",
        outdir + "/duplications.bedpe"
    message:
        "Aggregating the detected variants"
    run:
        from pathlib import Path
        with open(output[0], "w") as inversions, open(output[1], "w") as deletions, open(output[2], "w") as duplications:
            header = ["Sample","Chr1","Break1","Chr2","Break2","SplitMolecules","DiscordantReads","Orientation","Haplotype","Score","PassFilter","SV"]
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

rule genome_setup:
    input:
        genomefile
    output: 
        f"Genome/{bn}"
    container:
        None
    message: 
        "Symlinking {input}"
    shell: 
        """
        if (file {input} | grep -q compressed ) ;then
            # is regular gzipped, needs to be BGzipped
            seqtk seq {input} | bgzip -c > {output}
        else
            # if BZgipped or isn't compressed, just copied
            cp -f {input} {output}
        fi
        """

rule genome_faidx:
    input: 
        f"Genome/{bn}"
    output: 
        fai = f"Genome/{bn}.fai",
        gzi = f"Genome/{bn}.gzi" if genome_zip else []
    log:
        f"Genome/{bn}.faidx.log"
    params:
        genome_zip
    container:
        None
    message:
        "Indexing {input}"
    shell: 
        """
        if [ "{params}" = "True" ]; then
            samtools faidx --gzi-idx {output.gzi} --fai-idx {output.fai} {input} 2> {log}
        else
            samtools faidx --fai-idx {output.fai} {input} 2> {log}
        fi
        """

rule create_report:
    input:
        bedpe = outdir + "/bedpe/{sample}.bedpe",
        fai   = f"Genome/{bn}.fai"
    output:
        outdir + "/reports/{sample}.naibr.html"
    conda:
        f"{envdir}/r.yaml"
    message:
        "Creating report: {wildcards.sample}"
    script:
        "report/naibr.Rmd"

rule workflow_summary:
    default_target: True
    input:
        bedpe = collect(outdir + "/bedpe/{sample}.bedpe", sample = samplenames),
        bedpe_agg = collect(outdir + "/{sv}.bedpe", sv = ["inversions", "deletions","duplications"]),
        reports =  collect(outdir + "/reports/{sample}.naibr.html", sample = samplenames) if not skipreports else []
    message:
        "Summarizing the workflow: {output}"
    run:
        os.system(f"rm -rf {outdir}/naibrlog")
        argdict = process_args(extra)
        with open(outdir + "/workflow/sv.naibr.summary", "w") as f:
            _ = f.write("The harpy sv naibr workflow ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write("naibr variant calling ran using these configurations:\n")
            _ = f.write(f"    bam_file=BAMFILE\n")
            _ = f.write(f"    prefix=PREFIX\n")
            _ = f.write(f"    outdir=Variants/naibr/PREFIX\n")
            for i in argdict:
                _ = f.write(f"    {i}={argdict[i]}\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")