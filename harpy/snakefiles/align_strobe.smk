import os
import re

localrules: all
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
bn            = bn[:-3] if bn.lower().endswith(".gz") else bn
bn_r          = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
workflow_geno = f"workflow/reference/{bn}"
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
        geno = workflow_geno,
        fai = f"{workflow_geno}.fai"
    log:
        f"{workflow_geno}.preprocess.log"
    shell: 
        """
        {{
            seqtk seq {input} > {output.geno}
            samtools faidx --fai-idx {output.fai} {output.geno}
        }} 2> {log}
        """

rule align:
    input:
        workflow_geno,
        get_fq,
    output:  
        bam = temp("strobealign/{sample}.strobe.bam"),
        tmp = temp(directory("strobealign/{sample}_tmp"))
    log:
        "logs/strobealign/{sample}.strobealign.log"
    params: 
        static = "-N 2 -C" if illumina_old else "-N 2",
        RGid = lambda wc: f"--rg-id={wc.get('sample')}",
        RGsm = lambda wc: f"--rg=SM:{wc.get('sample')}",
        extra = extra
    threads:
        12
    resources:
        tmpdir = lambda wc: f"strobealign/{wc.sample}_tmp"
    conda:
        "envs/align.yaml"
    container:
        f"docker://pdimens/harpy:align_{VERSION}"
    shell:
        """
        mkdir -p {resources.tmpdir}
        {{
            strobealign {params} -t {threads} {input} |
            samtools collate -T {resources.tmpdir} -O -u -l 0 -
        }} 2> {log} > {output.bam}
        """
