containerized: "docker://pdimens/harpy:latest"

import os
import random
import logging

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake_log"])
    logger.addHandler(logfile_handler)

envdir = os.path.join(os.getcwd(), "workflow", "envs")
genome = config["inputs"]["genome"]
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
        temp(f"{outprefix}.simseq.genome.fa"),
        f"{outprefix}.refseq2simseq.SNP.vcf" if snp else [],
        f"{outprefix}.refseq2simseq.INDEL.vcf" if indel else [],
        f"{outprefix}.refseq2simseq.map.txt"
    log:
        f"logs/{outprefix}.log"
    params:
        prefix = f"{outprefix}",
        parameters = variant_params
    conda:
        f"{envdir}/simulations.yaml"
    shell:
        "simuG -refseq {input.geno} -prefix {params.prefix} {params.parameters} > {log}"

rule rename_haploid:
    input:
        fasta = f"{outprefix}.simseq.genome.fa",
        snpvcf = f"{outprefix}.refseq2simseq.SNP.vcf" if snp else [],
        indelvcf = f"{outprefix}.refseq2simseq.INDEL.vcf" if indel else [],
        mapfile = f"{outprefix}.refseq2simseq.map.txt"
    output:
        fasta = f"{outprefix}.fasta.gz",
        snpvcf = f"{outprefix}.snp.vcf" if snp else [],
        indelvcf = f"{outprefix}.indel.vcf" if indel else [],
        mapfile = f"{outprefix}.map"
    run:
        shell(f"bgzip -c {input.fasta} > {output.fasta}")
        os.rename(input.mapfile, output.mapfile)
        if input.snpvcf:
            os.rename(input.snpvcf, output.snpvcf)
        if input.indelvcf:
            os.rename(input.indelvcf, output.indelvcf)

rule diploid_snps:
    input:
        f"{outprefix}.snp.vcf"
    output:
        f"haplotype_1/{outprefix}.hap1.snp.vcf",
        f"haplotype_2/{outprefix}.hap2.snp.vcf"
    params:
        het = heterozygosity
    run:
        os.makedirs("haplotype_1", exist_ok = True)
        os.makedirs("haplotype_2", exist_ok = True)
        rng = random.Random(randomseed) if randomseed else random.Random()
        with open(input[0], "r") as in_vcf, open(output[0], "w") as hap1, open(output[1], "w") as hap2:
            for line in in_vcf:
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
        f"{outprefix}.indel.vcf"
    output:
        f"haplotype_1/{outprefix}.hap1.indel.vcf",
        f"haplotype_2/{outprefix}.hap2.indel.vcf"

rule simulate_diploid:
    input:
        snp_hap = f"haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}.snp.vcf" if snp else [],
        indel_hap = f"haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}.indel.vcf" if indel else [],
        geno = genome
    output:
        temp(f"haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}.simseq.genome.fa"),
        temp(f"haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}.refseq2simseq.map.txt"),
        temp(f"haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}.refseq2simseq.INDEL.vcf") if indel else [],
        temp(f"haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}.refseq2simseq.SNP.vcf") if snp else []
    log:
        f"logs/{outprefix}.hap{{haplotype}}.log"
    params:
        prefix = f"haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}",
        snp = f"-snp_vcf haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}.snp.vcf" if snp else "",
        indel = f"-indel_vcf haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}.indel.vcf" if indel else ""
    conda:
        f"{envdir}/simulations.yaml"
    shell:
        "simuG -refseq {input.geno} -prefix {params.prefix} {params.snp} {params.indel} > {log}"

rule rename_diploid:
    input:
        fasta= f"haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}.simseq.genome.fa",
        mapfile = f"haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}.refseq2simseq.map.txt",
    output:
        fasta = f"haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}.fasta.gz",
        mapfile = f"haplotype_{{haplotype}}/{outprefix}.hap{{haplotype}}.map"
    container:
        None
    shell:
        """
        bgzip -c {input.fasta} > {output.fasta}
        cp {input.mapfile} {output.mapfile}
        """

rule workflow_summary:
    default_target: True
    input:
        f"{outprefix}.fasta.gz",
        collect(f"{outprefix}" + ".{var}.vcf", var = variants),
        collect("haplotype_{n}/" + outprefix + ".hap{n}.fasta.gz", n = [1,2]) if heterozygosity > 0 and not only_vcf else [],
        collect("haplotype_{n}/" + outprefix + ".hap{n}" + ".{var}.vcf", n = [1,2], var = variants) if heterozygosity > 0 else [],
        collect("haplotype_{n}/" + outprefix + ".hap{n}" + ".map", n = [1,2]) if heterozygosity > 0 else []
    params:
        prefix = f"{outprefix}",
        parameters = variant_params,
        snp = f"-snp_vcf haplotype_X/{outprefix}.snp.hapX.vcf" if snp else "",
        indel = f"-indel_vcf haplotype_X/{outprefix}.indel.hapX.vcf" if indel else ""
    run:
        summary = ["The harpy simulate snpindel workflow ran using these parameters:"]
        summary.append(f"The provided genome: {genome}")
        summary.append(f"Heterozygosity specified: {heterozygosity}")
        haploid = "Haploid variants were simulated using simuG:\n"    
        haploid += f"\tsimuG -refseq {genome} -prefix {params.prefix} {params.parameters}"
        summary.append(haploid)
        if heterozygosity > 0 and not only_vcf:
            diploid = "Diploid variants were simulated after splitting by the heterozygosity ratio:\n"
            diploid += f"\tsimuG -refseq {genome} -prefix hapX {params.snp} {params.indel}"
            summary.append(diploid)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['snakemake_command']}"
        summary.append(sm)
        with open(f"workflow/simulate.snpindel.summary", "w") as f:
            f.write("\n\n".join(summary))