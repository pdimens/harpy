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
genome = config["inputs"]["genome"]
envdir = os.getcwd() + "/.harpy_envs"
snp_vcf = config["snp"].get("vcf", None)
indel_vcf = config["indel"].get("vcf", None)
heterozygosity = float(config["heterozygosity"]["ratio"])
only_vcf = config["heterozygosity"]["only_vcf"]
outprefix = config["prefix"]
randomseed = config.get("random_seed", None)
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
    snp_count = config["snp"].get("count", None)
    indel_count =  config["indel"].get("count", None)
    variant_params = ""
    if snp_count:
        snp = True
        variant_params += f" -snp_count {snp_count}"
        snp_constraint = config["snp"].get("gene_constraints", None)
        variant_params += f" -coding_partition_for_snp_simulation {snp_constraint}" if snp_constraint else ""
        ratio = config["snp"].get("titv_ratio", None)
        variant_params += f" -titv_ratio {ratio}" if ratio else ""
    if indel_count:
        indel = True
        variant_params += f" -indel_count {indel_count}"
        ratio = config["indel"].get("indel_ratio", None)
        variant_params += f" -ins_del_ratio {ratio}" if ratio else ""
        size_alpha = config["indel"].get("size_alpha", None)
        variant_params += f" -indel_size_powerlaw_alpha {size_alpha}" if size_alpha else ""        
        size_constant = config["indel"].get("size_constant", None)
        variant_params += f" -indel_size_powerlaw_constant {size_constant}" if size_constant else ""        

    centromeres = config["inputs"].get("centromeres", None)
    variant_params += f" -centromere_gff {centromeres}" if centromeres else ""
    genes = config["inputs"].get("genes", None)
    variant_params += f" -gene_gff {genes}" if genes else ""
    exclude = config["inputs"].get("excluded_chromosomes", None)
    variant_params += f" -excluded_chr_list {exclude}" if exclude else ""
    variant_params += f" -seed {randomseed}" if randomseed else ""

variants = [i for i,j in zip(["snp", "indel"], [snp, indel]) if j]

if snp_vcf:
    rule convert_snp_vcf:
        input:
            snp_vcf
        output:
            snp_vcf_correct
        container:
            None
        shell:
            "bcftools view -Oz {input} > {output}"

if indel_vcf:
    rule convert_indel_vcf:
        input:
            indel_vcf
        output:
            indel_vcf_correct
        container:
            None
        shell:
            "bcftools view -Oz {input} > {output}"

rule simulate_haploid:
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
        f"{envdir}/simulations.yaml"
    shell:
        "perl {params.simuG} -refseq {input.geno} -prefix {params.prefix} {params.parameters} > {log}"

rule diploid_snps:
    input:
        f"{outdir}/{outprefix}.snp.vcf"
    output:
        f"{outdir}/diploid/{outprefix}.snp.hap1.vcf",
        f"{outdir}/diploid/{outprefix}.snp.hap2.vcf"
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

use rule diploid_snps as diploid_indels with:
    input:
        f"{outdir}/{outprefix}.indel.vcf"
    output:
        f"{outdir}/diploid/{outprefix}.indel.hap1.vcf",
        f"{outdir}/diploid/{outprefix}.indel.hap2.vcf"

rule simulate_diploid:
    input:
        snp_hap = f"{outdir}/diploid/{outprefix}.snp.hap{{haplotype}}.vcf" if snp else [],
        indel_hap = f"{outdir}/diploid/{outprefix}.indel.hap{{haplotype}}.vcf" if indel else [],
        geno = genome
    output:
        f"{outdir}/diploid/{outprefix}.hap{{haplotype}}.fasta",
        temp(f"{outdir}/diploid/{outprefix}.hap{{haplotype}}.indel.vcf") if indel else [],
        temp(f"{outdir}/diploid/{outprefix}.hap{{haplotype}}.snp.vcf") if snp else []
    log:
        f"{outdir}/logs/{outprefix}.hap{{haplotype}}.log"
    params:
        prefix = f"{outdir}/diploid/{outprefix}.hap{{haplotype}}",
        simuG = f"{outdir}/workflow/scripts/simuG.pl",
        snp = f"-snp_vcf {outdir}/diploid/{outprefix}.snp.hap{{haplotype}}.vcf" if snp else "",
        indel = f"-indel_vcf {outdir}/diploid/{outprefix}.indel.hap{{haplotype}}.vcf" if indel else ""
    conda:
        f"{envdir}/simulations.yaml"
    shell:
        "perl {params.simuG} -refseq {input.geno} -prefix {params.prefix} {params.snp} {params.indel} > {log}"

rule workflow_summary:
    default_target: True
    input:
        multiext(f"{outdir}/{outprefix}", ".bed", ".fasta"),
        collect(f"{outdir}/{outprefix}" + ".{var}.vcf", var = variants),
        collect(f"{outdir}/diploid/{outprefix}" + ".hap{n}.fasta", n = [1,2]) if heterozygosity > 0 and not only_vcf else [],
        collect(f"{outdir}/diploid/{outprefix}" + ".{var}.hap{n}.vcf", n = [1,2], var = variants) if heterozygosity > 0 else []
    params:
        prefix = f"{outdir}/{outprefix}",
        parameters = variant_params,
        dip_prefix = f"{outdir}/diploid/{outprefix}.hap{{haplotype}}",
        snp = f"-snp_vcf {outdir}/diploid/{outprefix}.snp.hap{{haplotype}}.vcf" if snp else "",
        indel = f"-indel_vcf {outdir}/diploid/{outprefix}.indel.hap{{haplotype}}.vcf" if indel else ""
    run:
        summary = ["The harpy simulate snpindel workflow ran using these parameters:"]
        summary.append(f"The provided genome: {genome}")
        summary.append(f"Heterozygosity specified: {heterozygosity}")
        haploid = "Haploid variants were simulated using simuG:\n"    
        haploid += f"simuG -refseq {genome} -prefix {params.prefix} {params.parameters}"
        summary.append(haploid)
        if heterozygosity > 0 and not only_vcf:
            diploid = "Diploid variants were simulated after splitting by the heterozygosity ratio:\n"
            diploid += f"\tsimuG -refseq {genome} -prefix {params.dip_prefix} {params.snp} {params.indel}"
            summary.append(diploid)
        sm = "The Snakemake workflow was called via command line:"
        sm += f"\t{config["workflow_call"]}"
        summary.append(sm)
        with open(f"{outdir}/workflow/simulate.snpindel.summary", "w") as f:
            f.write("\n\n".join(summary))