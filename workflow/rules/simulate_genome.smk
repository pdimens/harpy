import sys
import random
from rich.panel import Panel
from rich import print as rprint

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

rule simulate_snps:
    input:
        genome
    output:
        expand(prefix.{ext}, ext = "SNP.vcf", "variants.bed", "simulated.fasta")
    params:
    
    conda:

    message:
        "Simulating inversions"
    shell:
        """
        perl simuG.pl -refseq ../dmel.trunc.fa -inversion_count 60 -inversion_min_size 1000 -inversion_max_size 25000000 -seed 6969696969 -prefix inversions
        """

rule simulate_snps:
    input:
        genome
    output:
        expand(prefix.{ext}, ext = "SNP.vcf", "variants.bed", "simulated.fasta")
    params:
    
    conda:

    message:
        "Simulating SNPs and indels"
    shell:
        """
        perl simuG.pl -refseq inversions.simseq.genome.fa -indel_count 100 -snp_count 180000 -translocation_count 0 -seed 6969696969 -prefix snp_inversions
        """

rule create_heterozygous_snp:
    input:
        "prefix.SNP.vcf"
    output:
        "prefix.SNP.hap2.vcf"
    params:
        het_threshold = 0.4
    message:
        "Creating VCF file of heterozygous SNPs to simulate"
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


rule create_heterozygous_inv:
    input:
        "prefix.inversion.vcf"
    output:
        "prefix.ivnersion.hap2.vcf"
    params:
        het_threshold = 0.4
    message:
        "Creating VCF file of heterozygous inversions to simulate"
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

rule create_heterozygous_cnv:
    input:
        "prefix.CNV.vcf"
    output:
        "prefix.CNV.hap2.vcf"
    params:
        het_threshold = 0.4
    message:
        "Creating VCF file of heterozygous CNVs to simulate"
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

rule create_heterozygous_translocation:
    input:
        "prefix.translocation.vcf"
    output:
        "prefix.translocation.hap2.vcf"
    params:
        het_threshold = 0.4
    message:
        "Creating VCF file of heterozygous translocations to simulate"
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