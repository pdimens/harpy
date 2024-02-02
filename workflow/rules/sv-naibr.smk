from rich import print as rprint
from rich.panel import Panel
import sys
import os
import re

bam_dir     = config["seq_directory"]
samplenames = config["samplenames"] 
extra       = config.get("extra", "") 
genomefile  = config["genomefile"]
molecule_distance = config["molecule_distance"]
outdir      = "Variants/naibr"
skipreports = config["skipreports"]
bn          = os.path.basename(genomefile)
genome_zip  = True if (bn.endswith(".gz") or bn.endswith(".GZ")) else False
bn_idx      = f"{bn}.gzi" if genome_zip else f"{bn}.fai"

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

def process_args(args):
    argsDict = {
        "min_mapq" : 30,
        "d"        : molecule_distance,
        "min_sv"   : 1000,
        "k"        : 3
    }
    if args != "":
        words = [i for i in re.split("\s|=", args) if len(i) > 0]
        for i in zip(words[::2], words[1::2]):
            argsDict[i[0]] = i[1]
    return argsDict

onerror:
    print("")
    rprint(
        Panel(
            f"The workflow has terminated due to an error. See the log file below for more details.",
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

rule index_alignment:
    input:
        bam_dir + "/{sample}.bam"
    output:
        bam_dir + "/{sample}.bam.bai"
    message:
        "Indexing alignment: {wildcards.sample}"
    shell:
        "samtools index {input} {output} 2> /dev/null"

rule create_config:
    input:
        bam_dir + "/{sample}.bam"
    output:
        outdir + "/configs/{sample}.config"
    params:
        lambda wc: wc.get("sample")
    message:
        "Creating naibr config file: {wildcards.sample}"
    run:
        argdict = process_args(extra)
        with open(output[0], "w") as conf:
            _ = conf.write(f"bam_file={input[0]}\n")
            _ = conf.write(f"prefix={params[0]}\n")
            _ = conf.write(f"outdir=Variants/naibr/{params[0]}\n")
            _ = conf.write(f"threads={workflow.cores}\n")
            for i in argdict:
                _ = conf.write(f"{i}={argdict[i]}\n")

rule call_sv:
    input:
        bam   = bam_dir + "/{sample}.bam",
        bai   = bam_dir + "/{sample}.bam.bai",
        conf  = outdir + "/configs/{sample}.config"
    output:
        bedpe = outdir + "/{sample}/{sample}.bedpe",
        refmt = outdir + "/{sample}/{sample}.reformat.bedpe",
        vcf   = outdir + "/{sample}/{sample}.vcf"
    log:
        full  = outdir + "/logs/{sample}.naibrlog",
        pairs = outdir + "/logs/{sample}.log"
    threads:
        8
    conda:
        os.getcwd() + "/harpyenvs/variants.sv.yaml"     
    message:
        "Calling variants: {wildcards.sample}"
    shell:
        """
        naibr {input.conf} > {log.full} 2>&1
        grep -v "pairs/s" {log.full} > {log.pairs}
        """

rule infer_sv:
    input:
        bedpe = outdir + "/{sample}/{sample}.bedpe",
        refmt = outdir + "/{sample}/{sample}.reformat.bedpe",
        vcf   = outdir + "/{sample}/{sample}.vcf"
    output:
        bedpe = outdir + "/{sample}.bedpe",
        refmt = outdir + "/IGV/{sample}.reformat.bedpe",
        fail  = outdir + "/filtered/{sample}.fail.bedpe",
        vcf   = outdir + "/vcf/{sample}.vcf" 
    params:
        outdir = lambda wc: outdir + "/" + wc.get("sample")
    message:
        "Inferring variants from naibr output: {wildcards.sample}"
    shell:
        """
        inferSV.py {input.bedpe} -f {output.fail} > {output.bedpe}
        mv {input.refmt} {output.refmt} &&
        mv {input.vcf} {output.vcf} &&
        rm -rf {params.outdir}
        """

rule genome_link:
    input:
        genomefile
    output: 
        f"Genome/{bn}"
    message: 
        "Symlinking {input}"
    shell: 
        """
        if (file {input} | grep -q compressed ) ;then
            # is regular gzipped, needs to be BGzipped
            zcat {input} | bgzip -c > {output}
        elif (file {input} | grep -q BGZF ); then
            # is bgzipped, just linked
            ln -sr {input} {output}
        else
            # isn't compressed, just linked
            ln -sr {input} {output}
        fi
        """

if genome_zip:
    rule genome_compressed_faidx:
        input: 
            f"Genome/{bn}"
        output: 
            gzi = f"Genome/{bn}.gzi",
            fai = f"Genome/{bn}.fai"
        log:
            f"Genome/{bn}.faidx.gzi.log"
        message:
            "Indexing {input}"
        shell: 
            "samtools faidx --gzi-idx {output.gzi} --fai-idx {output.fai} {input} 2> {log}"
else:
    rule genome_faidx:
        input: 
            f"Genome/{bn}"
        output: 
            f"Genome/{bn}.fai"
        log:
            f"Genome/{bn}.faidx.log"
        message:
            "Indexing {input}"
        shell:
            "samtools faidx --fai-idx {output} {input} 2> {log}"

rule report:
    input:
        bedpe = outdir + "/{sample}.bedpe",
        fai   = f"Genome/{bn}.fai"
    output:
        outdir + "/reports/{sample}.naibr.html"
    conda:
        os.getcwd() + "/harpyenvs/r-env.yaml"
    message:
        "Creating report: {wildcards.sample}"
    script:
        "report/Naibr.Rmd"

rule log_runtime:
    output:
        outdir + "/workflow/sv.naibr.workflow.summary"
    message:
        "Creating record of relevant runtime parameters: {output}"
    run:
        argdict = process_args(extra)
        with open(output[0], "w") as f:
            _ = f.write("The harpy variants sv module ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write(f"The directory with alignments: {bam_dir}\n\n")
            _ = f.write("naibr variant calling ran using these configurations:\n")
            _ = f.write(f"    bam_file=BAMFILE\n")
            _ = f.write(f"    prefix=PREFIX\n")
            _ = f.write(f"    outdir=Variants/naibr/PREFIX\n")
            for i in argdict:
                _ = f.write(f"    {i}={argdict[i]}\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]))

results = list()
results.apppend(expand(outdir + "/{sample}.bedpe", sample = samplenames))
results.apppend(outdir + "/workflow/sv.naibr.workflow.summary")
if not skipreports:
    results.apppend(expand(outdir + "/reports/{sample}.naibr.html", sample = samplenames))

rule all:
    default_target: True
    input:
    message:
        "Checking for expected workflow output"
    shell:
        "rm -rf Variants/naibrlog"
