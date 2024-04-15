import sys
import os
import random
from rich.panel import Panel
from rich import print as rprint

indir = config["input_directory"]
outdir = config["output_directory"]
genome = config["genome"]
snp_vcf = config.get("snp_vcf", None)
indel_vcf = config.get("indel_vcf", None)
heterozygosity = config["heterozygosity"]
outprefix = config["prefix"]
in_vcfs = []
snp = False 
indel = False

if snp_vcf or indel_vcf:
    variant_params = ""
    if snp_vcf:
        snp = True
        snp_vcf_correct = snp_vcf[:-4] + ".vcf.gz" if snp_vcf.lower().endswith("bcf") else snp_vcf
        variant_params += f" -snp_vcf {snp_vcf_correct}"
        in_vcfs.append(snp_vcf_correct)
    if indel_vcf:
        indel = True
        indel_vcf_correct = indel_vcf[:-4] + ".vcf.gz" if indel_vcf.lower().endswith("bcf") else indel_vcf
        variant_params += f" -indel_vcf {indel_vcf_correct}"
        in_vcfs.append(indel_vcf_correct)
else:
    snp_count = config.get("snp_count", None)
    indel_count =  config.get("indel_count", None)
    variant_params = ""
    if snp_count:
        snp = True
        variant_params += f" -snp_count {snp_count}"
        snp_constraint = config.get("snp_gene_constraints", None)
        variant_params += f" -coding_partition_for_snp_simulation {snp_constraint}" if snp_constraint else ""
        ratio = config.get("titv_ratio", None)
        variant_params += f" -titv_ratio {ratio}" if ratio else ""

    if indel_count:
        indel = True
        variant_params += f" -indel_count {indel_count}"
        ratio = config.get("indel_ratio", None)
        variant_params += f" -ins_del_ratio {ratio}" if ratio else ""

    centromeres = config.get("centromeres", None)
    variant_params += f" -centromere_gff {indir + '/' + os.path.basename(centromeres)}" if centromeres else ""
    genes = config.get("genes", None)
    variant_params += f" -gene_gff {indir + '/' + os.path.basename(genes)}" if genes else ""
    exclude = config.get("exclude_chr", None)
    variant_params += f" -excluded_chr_list {indir + '/' + os.path.basename(exclude)}" if exclude else ""
    randomseed = config.get("randomseed", None)
    variant_params += f" -seed {randomseed}" if randomseed else ""

variants = [i for i,j in zip(["snp", "indel"], [snp, indel]) if j]

onsuccess:
    print("")
    rprint(
        Panel(
            f"The workflow has finished successfully! Find the results in [bold]{outdir}/[/bold]",
            title = f"[bold]harpy simulate snpindel",
            title_align = "left",
            border_style = "green"
            ),
        file = sys.stderr
    )

onerror:
    print("")
    rprint(
        Panel(
            f"The workflow has terminated due to an error. See the log file below for more details.",
            title = "[bold]harpy simulate snpindel",
            title_align = "left",
            border_style = "red"
            ),
        file = sys.stderr
    )

if snp_vcf:
    rule convert_snpvcf:
        input:
            snp_vcf
        output:
            snp_vcf_correct
        message:
            "Converting {input} to compressed VCF format"
        shell:
            "bcftools view -Oz {input} > {output}"

if indel_vcf:
    rule convert_indelvcf:
        input:
            indel_vcf
        output:
            indel_vcf_correct
        message:
            "Converting {input} to compressed VCF format"
        shell:
            "bcftools view -Oz {input} > {output}"

rule simulate_variants:
    input:
        in_vcfs,
        geno = genome
    output:
        multiext(f"{outdir}/{outprefix}", ".bed", ".fasta"),
        collect(f"{outdir}/{outprefix}." + "{var}.vcf", var = variants)
    log:
        f"{outdir}/logs/{outprefix}.log"
    params:
        prefix = f"{outdir}/{outprefix}",
        simuG = f"{outdir}/workflow/scripts/simuG.pl",
        parameters = variant_params
    conda:
        os.getcwd() + "/.harpy_envs/simulations.yaml"
    message:
        "Simulating snps and/or indels"
    shell:
        """
        perl {params.simuG} -refseq {input.geno} -prefix {params.prefix} {params.parameters} > {log}
        """

rule create_heterozygote_snp_vcf:
    input:
        f"{outdir}/{outprefix}.snp.vcf"
    output:
        f"{outdir}/{outprefix}.snp.hap1.vcf",
        f"{outdir}/{outprefix}.snp.hap2.vcf"
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

rule create_heterozygote_indel_vcf:
    input:
        f"{outdir}/{outprefix}.indel.vcf"
    output:
        f"{outdir}/{outprefix}.indel.hap1.vcf",
        f"{outdir}/{outprefix}.indel.hap2.vcf"
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
        multiext(f"{outdir}/{outprefix}", ".bed", ".fasta"),
        collect(f"{outdir}/{outprefix}." + "{var}.vcf", var = variants),
        collect(f"{outdir}/{outprefix}" + ".{var}.hap{n}.vcf", n = [1,2], var = variants) if heterozygosity > 0 else []
    message:
        "Checking for workflow outputs"