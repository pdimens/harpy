containerized: "docker://pdimens/harpy:latest"

import os
import random
import logging

onstart:
    logger.logger.addHandler(logging.FileHandler(config["snakemake_log"]))
onsuccess:
    os.remove(logger.logfile)
onerror:
    os.remove(logger.logfile)

outdir = config["output_directory"]
envdir = os.getcwd() + "/.harpy_envs"
variant = config["variant_type"]
outprefix = config["prefix"]
genome = config["inputs"]["genome"]
vcf = config["inputs"].get("vcf", None)
heterozygosity = float(config["heterozygosity"]["value"])
only_vcf = config["heterozygosity"]["only_vcf"]
randomseed = config.get("random_seed", None)

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
        collect(f"{outdir}/{outprefix}" + "{ext}", ext = [".vcf", ".bed", ".fasta"])
    log:
        f"{outdir}/logs/{outprefix}.log"
    params:
        prefix = f"{outdir}/{outprefix}",
        simuG = f"{outdir}/workflow/scripts/simuG.pl",
        parameters = variant_params
    conda:
        f"{envdir}/simulations.yaml"
    shell:
        "perl {params.simuG} -refseq {input.geno} -prefix {params.prefix} {params.parameters} > {log}"

rule diploid_variants:
    input:
        f"{outdir}/{outprefix}.vcf"
    output:
        f"{outdir}/diploid/{outprefix}.{variant}.hap1.vcf",
        f"{outdir}/diploid/{outprefix}.{variant}.hap2.vcf"
    params:
        het = heterozygosity
    run:
        rng = random.Random(randomseed) if randomseed else random.Random()
        with open(input[0], "r") as in_vcf, open(output[0], "w") as hap1, open(output[1], "w") as hap2:
            while True:
                line = in_vcf.readline()
                if not line:
                    break
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
        hap = f"{outdir}/diploid/{outprefix}.{variant}.hap{{haplotype}}.vcf",
        geno = genome
    output:
        f"{outdir}/diploid/{outprefix}.hap{{haplotype}}.fasta",
        temp(f"{outdir}/diploid/{outprefix}.hap{{haplotype}}.vcf")
    log:
        f"{outdir}/logs/{outprefix}.hap{{haplotype}}.log"
    params:
        prefix = f"{outdir}/diploid/{outprefix}.hap{{haplotype}}",
        simuG = f"{outdir}/workflow/scripts/simuG.pl",
        vcf_arg = f"-{variant}_vcf"
    conda:
        f"{envdir}/simulations.yaml"
    shell:
        "perl {params.simuG} -refseq {input.geno} -prefix {params.prefix} {params.vcf_arg} {input.hap} > {log}"

rule workflow_summary:
    default_target: True
    input:
        multiext(f"{outdir}/{outprefix}", ".vcf", ".bed", ".fasta"),
        collect(f"{outdir}/diploid/{outprefix}.hap" + "{n}.fasta", n = [1,2]) if heterozygosity > 0 and not only_vcf else [],
        collect(f"{outdir}/diploid/{outprefix}.{variant}.hap" + "{n}.vcf", n = [1,2]) if heterozygosity > 0 else []
