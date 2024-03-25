import sys
import random
from rich.panel import Panel
from rich import print as rprint

haps = 1 if heterozygosity == 0 else [1,2]

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

rule simulate_variants:
    input:
        genome
    output:
        expand(f"{outdir}/simulation.hap1".{ext}, ext = f"{variant}.vcf", "variants.bed", "fasta")
    params:
        prefix = f"{outdir}/{variant}.hap1",
        simuG = f"{outdir}/workflow/scripts/simuG.pl",
        parameters = #TODO PARAMS HERE
    conda:
        os.getcwd() + "/.harpy_envs/simulations.yaml"
    message:
        f"Simulating {variant}s for first haplotype"
    shell:
        """
        perl {params.simuG} -refseq {input} -prefix {params.prefix} {params.parameters}
        """

rule create_heterozygous_snp:
    input:
        f"{prefix}.SNP.vcf"
    output:
        f"{prefix}.SNP.hap2.vcf"
    params:
        heterozygosity
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

rule simulate_variants_hap2:
    input:
        genome = f"{outdir}/simulation.hap1.fasta",
        vcf    = f"{outdir}/simulation.hap1.{variant}.vcf"
    output:
        expand(f"{outdir}/simulation.hap2.{ext}", ext = f"{variant}.vcf", "variants.bed", "fasta")
    params:
        prefix = f"{outdir}/{variant}.hap2",
        vcf_arg = variant.lower() + "_vcf",
        simuG = f"{outdir}/workflow/scripts/simuG.pl"
    conda:
        os.getcwd() + "/.harpy_envs/simulations.yaml"
    message:
        f"Simulating {variant}s for 2nd haplotype"
    shell:
        """
        perl {params.simuG} -refseq {input.genome} -prefix {params.prefix} -{params.vcf_arg} {input.vcf}
        """

rule all:
    input:
        expand(f"{outdir}/simulation.{hap}.{ext}", hap = haps, ext = f"{variant}.vcf", "variants.bed", "fasta")
    message:
        "Checking for workflow outputs"