containerized: "docker://pdimens/harpy:latest"

import os
import random
import logging

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake"]["log"])
    logger.addHandler(logfile_handler)

variant = config["workflow"].split("_")[-1]
simuG_variant = variant.upper() if variant == "cnv" else variant
outprefix = config["prefix"]
genome = config["inputs"]["genome"]
vcf = config[variant].get("vcf", None)
heterozygosity = float(config["heterozygosity"]["ratio"])
only_vcf = config["heterozygosity"]["only_vcf"]
randomseed = config.get("random_seed", None)

if vcf:
    vcf_correct = vcf[:-4] + ".vcf.gz" if vcf.lower().endswith("bcf") else vcf
    variant_params = f"-{variant}_vcf {vcf_correct}"
else:
    variant_params = f"-{variant}_count " + str(config[variant]["count"])
    centromeres = config["inputs"].get("centromeres", None)
    variant_params += f" -centromere_gff {centromeres}" if centromeres else ""
    genes = config["inputs"].get("genes", None)
    variant_params += f" -gene_gff {genes}" if genes else ""
    exclude = config["inputs"].get("excluded_chromosomes", None)
    variant_params += f" -excluded_chr_list {exclude}" if exclude else ""
    variant_params += f" -seed {randomseed}" if randomseed else ""
    if variant in ["inversion", "cnv"]:  
        variant_params += f" -{variant}_min_size " +  str(config[variant]["min_size"])
        variant_params += f" -{variant}_max_size " +  str(config[variant]["max_size"])
    if variant == "cnv":
        variant_params += f" -duplication_tandem_dispersed_ratio " +  str(config[variant]["duplication_ratio"])
        variant_params += f" --cnv_max_copy_number " +  str(config[variant]["max_copy"])
        variant_params += f" --cnv_gain_loss_ratio " +  str(config[variant]["gain_ratio"])

if vcf:
    rule convert_vcf:
        input:
            vcf
        output:
            vcf_correct
        container:
            None
        shell:
            "bcftools view -Oz {input} > {output}"

rule simulate_haploid:
    input:
        vcf_correct if vcf else [],
        geno = genome
    output:
        temp(f"{outprefix}.simseq.genome.fa"),
        f"{outprefix}.refseq2simseq.{simuG_variant}.vcf",
        f"{outprefix}.refseq2simseq.map.txt"
    log:
        f"logs/{outprefix}.log"
    params:
        prefix = f"{outprefix}",
        parameters = variant_params
    conda:
        "envs/simulations.yaml"
    shell:
        "simuG -refseq {input.geno} -prefix {params.prefix} {params.parameters} > {log}"

rule rename_haploid:
    input:
        fasta = f"{outprefix}.simseq.genome.fa",
        vcf = f"{outprefix}.refseq2simseq.{simuG_variant}.vcf",
        mapfile = f"{outprefix}.refseq2simseq.map.txt"
    output:
        fasta = f"{outprefix}.fasta.gz",
        vcf = f"{outprefix}.{variant}.vcf",
        mapfile = f"{outprefix}.{variant}.map"
    run:
        shell(f"bgzip -c {input.fasta} > {output.fasta}")
        os.rename(input.mapfile, output.mapfile)
        os.rename(input.vcf, output.vcf)

rule diploid_variants:
    input:
        vcf = f"{outprefix}.{variant}.vcf"
    output:
        hap1 = f"haplotype_1/{outprefix}.hap1.{variant}.vcf",
        hap2 = f"haplotype_2/{outprefix}.hap2.{variant}.vcf"
    params:
        het = heterozygosity
    run:
        os.makedirs("haplotype_1", exist_ok = True)
        os.makedirs("haplotype_2", exist_ok = True)
        rng = random.Random(randomseed) if randomseed else random.Random()
        with open(input.vcf, "r") as in_vcf, open(output.hap1, "w") as hap1, open(output.hap2, "w") as hap2:
            for line in in_vcf:
                if line.startswith("#") or rng.uniform(0, 1) >= params.het:
                    # write header lines and homozygous variants to both files
                    hap1.write(line)
                    hap2.write(line)
                elif rng.random() < 0.5:
                    hap1.write(line)
                else:
                    hap2.write(line)

rule simulate_diploid:
    input:
        hap = f"haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}.{variant}.vcf",
        geno = genome
    output:
        temp(f"haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}.simseq.genome.fa"),
        temp(f"haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}.refseq2simseq.map.txt"),
        temp(f"haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}.refseq2simseq.{simuG_variant}.vcf")
    log:
        f"logs/{outprefix}.hap{{haplotype}}.log"
    params:
        prefix = f"haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}",
        vcf_arg = f"-{variant}_vcf"
    conda:
        "envs/simulations.yaml"
    shell:
        "simuG -refseq {input.geno} -prefix {params.prefix} {params.vcf_arg} {input.hap} > {log}"

rule rename_diploid:
    input:
        fasta = f"haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}.simseq.genome.fa",
        mapfile = f"haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}.refseq2simseq.map.txt"
    output:
        fasta = f"haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}.fasta.gz",
        mapfile = f"haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}.{variant}.map"
    container:
        None
    shell:
        """
        bgzip -c {input.fasta} > {output.fasta}
        cp {input.mapfile} {output.mapfile}
        """

rule all:
    default_target: True
    input:
        f"{outprefix}.fasta.gz",
        f"{outprefix}.{variant}.vcf",
        collect("haplotype_{n}/" + outprefix + ".hap{n}.fasta.gz", n = [1,2]) if heterozygosity > 0 and not only_vcf else [],
        collect("haplotype_{n}/" + outprefix + ".hap{n}" + f".{variant}.vcf", n = [1,2]) if heterozygosity > 0 else [],
        collect("haplotype_{n}/" + outprefix + ".hap{n}" + f".{variant}.map", n = [1,2]) if heterozygosity > 0 else []
