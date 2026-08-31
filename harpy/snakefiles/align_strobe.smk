import os
import re

wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

WORKFLOW   = config.get('Workflow') or {}
PARAMETERS = config.get('Parameters') or {}
REPORTS    = WORKFLOW.get("reports") or {} 
INPUTS     = config['Inputs']
VERSION    = WORKFLOW.get('harpy-version', 'latest')

illumina_old      = PARAMETERS.get("illumina-format-old", False)
extra 		      = PARAMETERS.get("extra", "") 
fqlist            = INPUTS["fastq"]
genomefile 	      = INPUTS["reference"]

bn 			  = os.path.basename(genomefile)
bn_r          = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
samplenames   = {re.sub(bn_r, "", os.path.basename(i), flags = re.IGNORECASE) for i in fqlist}
d             = dict(zip(samplenames, samplenames))

def get_fq(wildcards):
    # returns a list of fastq files for reads 1 and 2 based on *wildcards.sample* e.g.
    r = re.compile(fr".*/({re.escape(wildcards.sample)}){bn_r}", flags = re.IGNORECASE)
    return sorted(list(filter(r.match, fqlist))[:2])

rule process_reference:
    input:
        genomefile
    output: 
        geno = "workflow/reference/ref.fa.gz",
        fai = "workflow/reference/ref.fa.gz.fai"
    log:
        f"{bn}.preprocess.log"
    shell: 
        """
        {{
            seqtk seq {input} > {output.geno}
            samtools faidx --fai-idx {output.fai} {output.geno}
        }} 2> {log}
        """

rule align:
    input:
        "workflow/reference/ref.fa.gz",
        get_fq,
    output:  
        bam = temp("strobe/{sample}.strobe.bam"),
        tmp = temp(directory("strobe/{sample}_tmp"))
    log:
        "logs/strobealign/{sample}.strobe.log"
    params: 
        static = "-N 2 -C" if illumina_old else "-N 2",
        RGid = lambda wc: f"--rg-id={wc.get('sample')}",
        RGsm = lambda wc: f"--rg=SM:{wc.get('sample')}",
        extra = extra
    threads:
        12
    resources:
        tmpdir = lambda wc: f"strobe/{wc.sample}_tmp"
    conda:
        "envs/align.yaml"
    container:
        f"docker://pdimens/harpy:align_{VERSION}"
    shell:
        """
        mkdir -p {resources.tmpdir}
        {{
            strobealign {params} -t {threads} {input} |
            samtools collate -T {resources.tmpdir}/{wildcards.sample} -O -u -l 0 -
        }} 2> {log} > {output.bam}
        """
