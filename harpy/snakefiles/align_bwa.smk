import os
import re

wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

WORKFLOW    = config.get('Workflow') or {}
PARAMETERS = config.get('Parameters') or {}
REPORTS    = WORKFLOW.get("reports") or {} 
INPUTS     = config['Inputs']
VERSION    = WORKFLOW.get('harpy-version', 'latest')

illumina_old      = PARAMETERS.get("illumina-format-old", False)
extra 		      = PARAMETERS.get("extra", "") 
windowsize        = PARAMETERS.get("depth-windowsize", 50000)
fqlist            = INPUTS["fastq"]
genomefile 	      = INPUTS["reference"]

bn 			  = os.path.basename(genomefile)
bn_r          = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
bn_re = re.compile(bn_r, flags=re.IGNORECASE)
fq_by_sample = {}
for f in fqlist:
    name = bn_re.sub("", os.path.basename(f), count=1)
    fq_by_sample.setdefault(name, []).append(f)

samplenames = set(fq_by_sample)

def get_fq(wildcards):
    return sorted(fq_by_sample[wildcards.sample])[:2]


rule process_reference:
    input:
        genomefile
    output: 
        geno = "workflow/reference/ref.fa.gz",
        bwa_idx = multiext("workflow/reference/ref.fa.gz", ".l2b", ".mbw"),
        fai = "workflow/reference/ref.fa.gz.fai",
        gzi = "workflow/reference/ref.fa.gz.gzi"
    log:
        f"{bn}.preprocess.log"
    threads:
        4
    conda:
        "envs/align.yaml"
    container:
        f"docker://pdimens/harpy:align_{VERSION}"
    shell: 
        """
        {{
            seqtk seq {input} | bgzip -c > {output.geno}
            samtools faidx --gzi-idx {output.gzi} --fai-idx {output.fai} {output.geno}
            minibwa index -t {threads} {output.geno} 
        }} 2> {log}
        """

rule align:
    input:
        multiext("workflow/reference/ref.fa.gz", ".l2b", ".mbw"),
        ref   = "workflow/reference/ref.fa.gz",
        fastq = get_fq
    output:
        temp("bwa/{sample}.bwa.bam")
    log:
        "logs/bwa/{sample}.bwa.log"
    params:
        RG_tag = lambda wc: "-R \"@RG\\tID:" + wc.get("sample") + "\\tSM:" + wc.get("sample") + "\"",
        static = "-x sr -y" if illumina_old else "-x sr",
        extra = extra
    threads:
        12
    resources:
        tmpdir = lambda wc: f"bwa/{wc.sample}_tmp"
    conda:
        "envs/align.yaml"
    container:
        f"docker://pdimens/harpy:align_{VERSION}"
    shell:
        """
        mkdir -p {resources.tmpdir}; trap 'rm -rf {resources.tmpdir}' 0
        {{
            minibwa map -t {threads} {params} {input.ref} {input.fastq} |
            samtools collate -T {resources.tmpdir}/{wildcards.sample} -O -u - 
        }} 2> {log} > {output}
        """
