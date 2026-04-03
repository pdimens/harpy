import os
import re
from pathlib import Path

localrules: all, log_phasing
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

WORKFLOW   = config.get('Workflow') or {}
PARAMETERS = config.get('Parameters') or {}
INPUTS     = config['Inputs']
VERSION    = WORKFLOW.get('harpy-version', 'latest')

genomefile  = INPUTS["reference"]
bamlist     = INPUTS["alignments"]
vcffile     = INPUTS["vcf"]
samplenames = {Path(i).stem for i in bamlist}
mol_dist    = PARAMETERS.get("distance-threshold", 100000)
extra       = PARAMETERS.get("extra", None) 
ploidy       = PARAMETERS.get("ploidy", 2) 

bn            = os.path.basename(genomefile)
workflow_geno = f"workflow/reference/{bn[:-3]}" if bn.lower().endswith(".gz") else f"workflow/reference/{bn}"
vcfindex      = f"{vcffile}.csi" if vcffile.lower().endswith("bcf") else f"{vcffile}.tbi"

def get_alignments(wildcards):
    """returns a list with the bam file for the sample based on wildcards.sample"""
    r = re.compile(fr".*/({wildcards.sample})\.(bam|sam)$", flags = re.IGNORECASE)
    aln = list(filter(r.match, bamlist))
    return aln[0]

rule process_reference:
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

rule filter_invalid:
    input:
        get_alignments
    output:
        temp("filtered/{sample}.bam.bai"),
        invalid = temp("filtered/{sample}.invalid.bam"),    
        valid = temp("filtered/{sample}.bam")
    log:
        "logs/{sample}.filter_invalid.log"
    threads:
        4
    shell:
        """
        {{
            djinn sam filter-invalid -t {threads} --invalid {output.invalid} {input} > {output.valid}
            samtools index {output.valid}  
        }} 2> {log}
        """

rule index_vcf:
    input:
        vcffile
    output:
        vcffile + ".csi"
    shell:
        "bcftools index {input}"

rule index_vcf_gz:
    input:
        vcffile
    output:
        vcffile + ".tbi"
    shell:
        "tabix {input}"

rule phase_alignments:
    input:
        "filtered/{sample}.bam.bai",
        vcfindex,
        f"{workflow_geno}.fai",
        vcf = vcffile,
        aln = "filtered/{sample}.bam",
        ref = workflow_geno
    output:
        bam = temp("phased/{sample}.phased.bam"),
        log = "logs/{sample}.phase.log"
    log:
        "logs/{sample}.phase.log"
    params:
        f"--ploidy {ploidy}",
        f"-d {mol_dist}",
        "--tag-supplementary copy-primary",
        "--no-supplementary-strand-match",
        f"--supplementary-distance {3 * mol_dist}",
        "--ignore-read-groups",
        "--skip-missing-contigs",
        extra
    conda:
        "envs/phase.yaml"
    container:
        f"docker://pdimens/harpy:phase_{VERSION}"
    threads:
        3
    shell:
        """
        {{
            whatshap haplotag --sample {wildcards.sample} --reference {input.ref} {input.vcf} {input.aln} |
            samtools view -@ 1 -h -O BAM
        }} > {output.bam} 2> {output.log}
        """

rule log_phasing:
    input:
        collect("logs/{sample}.phase.log", sample = samplenames)
    output:
        "logs/phasing.summary.log"
    shell:
        """
        {{
            echo -e "sample\\tvalid_alignments\\tphased_alignments"
            for i in {input}; do
                SAMP=$(basename $i .phase.log)
                echo -e "${{SAMP}}\\t$(grep "Total alignments" $i)\\t$(grep "could be tagged" $i)" |
                    sed 's/ \\+ /\\t/g' | cut -f1,3,5
            done
        }} > {output}
        """

rule restore_invalid:
    input:
        "phased/{sample}.phased.bam",
        "filtered/{sample}.invalid.bam"
    output:
        pipe("{sample}.phased.unsort.sam")
    log:
        "logs/{sample}.restore_invalid.log"
    shell:
        "samtools merge -O SAM {output} {input} 2> {log}"
        
rule sort_phased_bam:
    input:
        "{sample}.phased.unsort.sam"
    output:
        "{sample}.phased.bam"
    log:
        "logs/{sample}.sort.log"
    resources:
        mem_mb = 2000
    threads:
        2
    shell:
        "samtools sort -@ 1 -o {output} -O BAM -m {resources.mem_mb}M {input} 2> {log}"

rule all:
    default_target: True
    input:
        bedpe = collect("{sample}.phased.bam", sample = samplenames),
        phaselog = "logs/phasing.summary.log"
