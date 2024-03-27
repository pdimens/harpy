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
    variant_params = f"-{variant}_vcf"
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
    if variant == "snp":
        snp_constraint = config.get("snp_gene_constraints", None)
        variant_params += f" -coding_partition_for_snp_simulation {snp_constraint}" if snp_constraint else ""
        ratio = config.get("ratio", None)
        variant_params += f" -titv_ratio {ratio}" if ratio else ""
    elif variant_type == "indel":
        ratio = config.get("ratio", None)
        variant_params += f" -ins_del_ratio {ratio}" if ratio else ""
    elif variant_type == "inversion":
        minsize = config.get("min_size", None)
        variant_params += f" -{variant}_min_size {min_size}" if min_size else ""
        maxsize = config.get("max_size", None)
        variant_params += f" -{variant}_max_size {max_size}" if max_size else ""
    elif variant_type == "cnv":
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
            title = "[bold]harpy simulate genome",
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
            title = "[bold]harpy simulate genome",
            title_align = "left",
            border_style = "red"
            ),
        file = sys.stderr
    )

#TODO this needs more logic before it and check the BCFTOOLS call
rule convert_vcf:
    input:
        f"{indir}/{vcf}"
    output:
        f"{indir}/{vcf_correct}"
    message:
        "Converting {input} to compressed VCF format"
    shell:
        "bcftools view -Oz {input} > {output}"

if vcf:
    rule simulate_variants:
        input:
            geno = f"{indir}/{genome}",
            vcf = vcf_correct
        output:
            expand(f"{outdir}/simulation.hap1".{ext}, ext = f"{variant}.vcf", "variants.bed", "fasta")
        params:
            prefix = f"{outdir}/{variant}.hap1",
            simuG = f"{outdir}/workflow/scripts/simuG.pl",
            parameters = variant_params
        conda:
            os.getcwd() + "/.harpy_envs/simulations.yaml"
        message:
            f"Simulating {variant}s for first haplotype"
        shell:
            """
            perl {params.simuG} -refseq {input.geno} -prefix {params.prefix} {params.parameters} {input.vcf}
            """
else:
    rule simulate_variants:
        input:
            genome
        output:
            expand(f"{outdir}/simulation.hap1".{ext}, ext = f"{variant}.vcf", "variants.bed", "fasta")
        params:
            prefix = f"{outdir}/{variant}.hap1",
            simuG = f"{outdir}/workflow/scripts/simuG.pl",
            parameters = variant_params
        conda:
            os.getcwd() + "/.harpy_envs/simulations.yaml"
        message:
            f"Simulating {variant}s for first haplotype"
        shell:
            """
            perl {params.simuG} -refseq {input} -prefix {params.prefix} {params.parameters}
            """

rule create_heterozygote_vcf:
    input:
        f"{outdir}.{variant}.vcf"
    output:
        f"{outdir}.{variant}.hap2.vcf"
    params:
        het
    message:
        "Creating subsampled variant file"
    run:
        random.seed(6969)
        out_vcf = open(output[0], "w")
        with open(input[0], "r") as in_vcf:
            while True:
                line = in_vcf.readline()
                if not line:
                    break
                if line.startswith("#"):
                    print(line.rstrip("\n"), file = out_vcf)
                else:
                    if random.uniform(0, 1) >= params[0]:
                        print(line.rstrip("\n"), file = out_vcf)

rule all:
    input:
        expand(f"{outdir}/simulation.hap1.{ext}", ext = f"{variant}.vcf", "variants.bed", "fasta"),
        f"{prefix}.{variant}.hap2.vcf" if heterozygosity > 0 else []
    message:
        "Checking for workflow outputs"