import sys
import os
import random
from rich.panel import Panel
from rich import print as rprint

indir = config["input_directory"] +
outdir = config["output_directory"]
variant = config["variant_type"]
genome = config["genome"]
vcf = config.get("vcf", None)
het = config["heterozygosity"]
if vcf:
    vcf_correct = vcf[:-4] + ".vcf.gz" if vcf.lower().endswith("bcf") else vcf
    variant_params = f"-{variant}_vcf {indir}/{vcf_correct}"
"
else:
    variant_params = f"-{variant}_count " + str(config["count"])
    centromeres = config.get("centromeres", None)
    variant_params += f" -centromere_gff {indir + '/' + os.path.basename(centromeres)}" if centromeres else ""
    genes = config.get("genes", None)
    variant_params += f" -gene_gff {indir + '/' + os.path.basename(genes)}" if genes else ""
    exclude = config.get("exclude_chr", None)
    variant_params += f" -excluded_chr_list {indir + '/' + os.path.basename(exclude)}" if exclude else ""
    randomseedseed = config.get("randomseed", None)
    variant_params += f" -seed {randomseed}" if randomseed else ""
    if variant_type == "inversion":
        minsize = config.get("min_size", None)
        variant_params += f" -{variant}_min_size {min_size}" if min_size else ""
        maxsize = config.get("max_size", None)
        variant_params += f" -{variant}_max_size {max_size}" if max_size else ""
    if variant_type == "cnv":
        minsize = config.get("min_size", None)
        variant_params += f" -{variant}_min_size {min_size}" if min_size else ""
        maxsize = config.get("max_size", None)
        variant_params += f" -{variant}_max_size {max_size}" if max_size else ""
        ratio   = config.get("ratio", None)
        variant_params += f" -duplication_tandem_dispersed_ratio {ratio}" if ratio else ""
        cnv_copy = config.get("cnv_max_copy", None)
        variant_params += f" --cnv_max_copy_number {cnv_copy}" if cnv_copy else ""
        cnv_ratio =config.get("cnv_ratio", None)
        variant_params += f" --cnv_gain_loss_ratio {cnv_ratio}" if cnv_ratio else ""

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
            f"The workflow has terminated due to an error. See the log file below for more details.",
            title = f"[bold]harpy simulate {variant}",
            title_align = "left",
            border_style = "red"
            ),
        file = sys.stderr
    )

rule convert_vcf:
    input:
        f"{indir}/{vcf}"
    output:
        f"{indir}/{vcf_correct}"
    message:
        "Converting {input} to compressed VCF format"
    shell:
        "bcftools view -Oz {input} > {output}"

rule simulate_variants:
    input:
        geno = genome,
        vcf_correct if vcf else []
    output:
        expand(f"{outdir}/{outprefix}.{variant}" + "{ext}", ext = [".vcf", ".bed", ".fasta"])
    params:
        prefix = f"{outdir}/{outprefix}.{variant}",
        simuG = f"{outdir}/workflow/scripts/simuG.pl",
        parameters = variant_params
    conda:
        os.getcwd() + "/.harpy_envs/simulations.yaml"
    message:
        f"Simulating {variant}s for first haplotype"
    shell:
        """
        perl {params.simuG} -refseq {input.geno} -prefix {params.prefix} {params.parameters}
        """

rule create_heterozygote_vcf:
    input:
        f"{outdir}/{outprefix}.{variant}.vcf"
    output:
        f"{outdir}/{outprefix}.{variant}.hap1.vcf",
        f"{outdir}/{outprefix}.{variant}.hap2.vcf"
    params:
        het
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
    input:
        expand(f"{outdir}/{outprefix}.{variant}" + "{ext}", ext = [".vcf", ".bed", ".fasta"]),
        expand(f"{outdir}/{prefix}.{variant}" + ".hap{n}.vcf", n = [1,2]) if heterozygosity > 0 else []
    message:
        "Checking for workflow outputs"