import os
import re
import logging
from pathlib import Path

onstart:
    logfile_handler = logger_manager._default_filehandler(config["Workflow"]["snakemake"]["log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

genomefile  = config["Inputs"]["reference"]
bamlist     = config["Inputs"]["alignments"]
vcffile     = config["Inputs"]["vcf"]
samplenames = {Path(i).stem for i in bamlist}
extra       = config["Parameters"].get("extra", None) 
mol_dist    = config["Parameters"]["molecule-distance"]
bn          = os.path.basename(genomefile)
if bn.lower().endswith(".gz"):
    workflow_geno = f"workflow/reference/{bn[:-3]}"
else:
    workflow_geno = f"workflow/reference/{bn}"
if vcffile.lower().endswith("bcf"):
    vcfindex = vcffile + ".csi"
else:
    vcfindex = vcffile + ".tbi"

def get_alignments(wildcards):
    """returns a list with the bam file for the sample based on wildcards.sample"""
    r = re.compile(fr".*/({wildcards.sample})\.(bam|sam)$", flags = re.IGNORECASE)
    aln = list(filter(r.match, bamlist))
    return aln[0]

rule preprocess_reference:
    input:
        genomefile
    output: 
        geno = workflow_geno,
        fai = f"{workflow_geno}.fai"
    log:
        f"{workflow_geno}.preprocess.log"
    shell: 
        """
        {{
            seqtk seq {input} > {output.geno}
            samtools faidx --fai-idx {output.fai} {output.geno}
        }} > {log}
        """

rule index_alignments:
    input:
        get_alignments
    output:
        temp("workflow/input/bam/{sample}.bam.bai"),
        bam = temp("workflow/input/bam/{sample}.bam")
    shell:
        """
        ln -sr {input} {output.bam}
        samtools index {output.bam}
        """

rule index_snps:
    input:
        vcffile
    output:
        vcffile + ".csi"
    shell:
        "bcftools index {input}"

rule index_snps_gz:
    input:
        vcffile
    output:
        vcffile + ".tbi"
    shell:
        "tabix {input}"

rule phase_alignments:
    input:
        "workflow/input/bam/{sample}.bam.bai",
        vcfindex,
        f"{workflow_geno}.fai",
        vcf = vcffile,
        aln = "workflow/input/bam/{sample}.bam",
        ref = workflow_geno
    output:
        bam = "{sample}.phased.bam",
        log = "logs/{sample}.phase.log"
    params:
        f"--linked-read-distance-cutoff {moldist}",
        "--tag-supplementary copy-primary",
        "--no-supplementary-strand-match",
        f"--supplementary-distance {moldist}",
        "--ignore-read-groups",
        "--skip-missing-contigs"
    conda:
        "envs/phase.yaml"
    container:
        "docker://pdimens/harpy:phase_latest"
    threads:
        4
    shell:
        "whatshap haplotag --sample {wildcards.sample} --output-threads={threads} -o {output.bam} --reference {input.ref} {input.vcf} {input.aln} 2> {output.log}"

rule log_phasing:
    input:
        collect("logs/{sample}.phase.log", sample = samplenames)
    output:
        "logs/whatshap-haplotag.log"
    shell:
        """
        {{
            echo -e "sample\\ttotal_alignments\\tphased_alignments"
            for i in {input}; do
                SAMP=$(basename $i .phaselog)
                echo -e "${{SAMP}}\\t$(grep "Total alignments" $i)\\t$(grep "could be tagged" $i)" |
                    sed 's/ \\+ /\\t/g' | cut -f1,3,5
            done
        }} > {output}
        """

rule all:
    default_target: True
    input:
        bedpe = collect("{sample}.phased.bam", sample = samplenames),
        phaselog = "logs/whatshap-haplotag.log"
