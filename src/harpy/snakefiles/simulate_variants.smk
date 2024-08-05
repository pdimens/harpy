containerized: "docker://pdimens/harpy:latest"

import sys
import os
import random
import logging as pylogging
from datetime import datetime
from rich.panel import Panel
from rich import print as rprint

outdir = config["output_directory"]
envdir = os.getcwd() + "/.harpy_envs"
variant = config["variant_type"]
outprefix = config["prefix"]
genome = config["inputs"]["genome"]
vcf = config["inputs"].get("vcf", None)
heterozygosity = float(config["heterozygosity"])
vcf_correct = "None"
if vcf:
    vcf_correct = vcf[:-4] + ".vcf.gz" if vcf.lower().endswith("bcf") else vcf
    variant_params = f"-{variant}_vcf {vcf_correct}"

else:
    variant_params = f"-{variant}_count " + str(config["count"])
    centromeres = config.get("centromeres", None)
    variant_params += f" -centromere_gff {centromeres}" if centromeres else ""
    genes = config["inputs"].get("genes", None)
    variant_params += f" -gene_gff {genes}" if genes else ""
    exclude = config["inputs"].get("exclude_chr", None)
    variant_params += f" -excluded_chr_list {exclude}" if exclude else ""
    randomseed = config.get("randomseed", None)
    variant_params += f" -seed {randomseed}" if randomseed else ""
    if variant == "inversion":
        min_size = config.get("min_size", None)
        variant_params += f" -{variant}_min_size {min_size}" if min_size else ""
        max_size = config.get("max_size", None)
        variant_params += f" -{variant}_max_size {max_size}" if max_size else ""
    if variant == "cnv":
        min_size = config.get("min_size", None)
        variant_params += f" -{variant}_min_size {min_size}" if min_size else ""
        max_size = config.get("max_size", None)
        variant_params += f" -{variant}_max_size {max_size}" if max_size else ""
        ratio   = config.get("ratio", None)
        variant_params += f" -duplication_tandem_dispersed_ratio {ratio}" if ratio else ""
        cnv_copy = config.get("cnv_max_copy", None)
        variant_params += f" --cnv_max_copy_number {cnv_copy}" if cnv_copy else ""
        cnv_ratio =config.get("cnv_ratio", None)
        variant_params += f" --cnv_gain_loss_ratio {cnv_ratio}" if cnv_ratio else ""

onstart:
    os.makedirs(f"{outdir}/logs/snakemake", exist_ok = True)
    dt_string = datetime.now().strftime("%d_%m_%Y-%H_%M_%S")
    extra_logfile_handler = pylogging.FileHandler(f"{outdir}/logs/snakemake/{dt_string}.snakelog")
    logger.logger.addHandler(extra_logfile_handler)

onsuccess:
    print("")
    rprint(
        Panel(
            f"The workflow has finished successfully! Find the results in [bold]{outdir}/[/bold]",
            title = f"[bold]harpy simulate {variant}",
            title_align = "left",
            border_style = "green"
            ),
        file = sys.stderr
    )

onerror:
    print("")
    rprint(
        Panel(
            f"The workflow has terminated due to an error. See the log file for more details:\n[bold]{outdir}/logs/snakemake/{dt_string}.snakelog[/bold]",
            title = f"[bold]harpy simulate {variant}",
            title_align = "left",
            border_style = "red"
            ),
        file = sys.stderr
    )

if vcf:
    rule convert_vcf:
        input:
            vcf
        output:
            vcf_correct
        container:
            None
        message:
            "Converting {input} to compressed VCF format"
        shell:
            "bcftools view -Oz {input} > {output}"

rule simulate_variants:
    input:
        vcf_correct if vcf else [],
        geno = genome
    output:
        collect(f"{outdir}/{outprefix}" + "{ext}", ext = [".vcf", ".bed", ".fasta"])
    log:
        f"{outdir}/logs/{outprefix}.log"
    params:
        prefix = f"{outdir}/{outprefix}",
        simuG = f"{outdir}/workflow/scripts/simuG.pl",
        parameters = variant_params
    conda:
        f"{envdir}/simulations.yaml"
    message:
        f"Simulating {variant}s onto genome"
    shell:
        "perl {params.simuG} -refseq {input.geno} -prefix {params.prefix} {params.parameters} > {log}"

rule create_heterozygote_vcf:
    input:
        f"{outdir}/{outprefix}.vcf"
    output:
        f"{outdir}/{outprefix}.hap1.vcf",
        f"{outdir}/{outprefix}.hap2.vcf"
    params:
        heterozygosity
    message:
        "Creating diploid variant files for heterozygosity = {params}"
    run:
        random.seed(6969)
        hap1_vcf, hap2_vcf = open(output[0], "w"), open(output[1], "w")
        with open(input[0], "r") as in_vcf:
            while True:
                line = in_vcf.readline()
                if not line:
                    break
                if line.startswith("#"):
                    hap1_vcf.write(line)
                    hap2_vcf.write(line)
                    continue
                if random.uniform(0, 1) >= params[0]:
                    # write homozygote
                    hap1_vcf.write(line)
                    hap2_vcf.write(line)
                else:
                    # 50% chance of falling into hap1 or hap2
                    if random.uniform(0, 1) >= 0.5:
                        hap1_vcf.write(line)
                    else:
                        hap2_vcf.write(line)

rule all:
    default_target: True
    input:
        multiext(f"{outdir}/{outprefix}", ".vcf", ".bed", ".fasta"),
        collect(f"{outdir}/{outprefix}.hap" + "{n}.vcf", n = [1,2]) if heterozygosity > 0 else []
    message:
        "Checking for workflow outputs"