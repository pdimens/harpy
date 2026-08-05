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
workflow_geno = f"workflow/reference/{bn}"
genome_zip    = True if bn.lower().endswith(".gz") else False
geno_idx      = f"{workflow_geno}.gzi" if genome_zip else f"{workflow_geno}.fai"
bn_r          = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
samplenames   = {re.sub(bn_r, "", os.path.basename(i), flags = re.IGNORECASE) for i in fqlist}
d             = dict(zip(samplenames, samplenames))

def get_fq(wildcards):
    # returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    r = re.compile(fr".*/({re.escape(wildcards.sample)}){bn_r}", flags = re.IGNORECASE)
    return sorted(list(filter(r.match, fqlist))[:2])

rule process_reference:
    input:
        genomefile
    output: 
        geno = workflow_geno,
        bwa_idx = multiext(workflow_geno, ".l2b", ".mbw"),
        fai = f"{workflow_geno}.fai",
        gzi = f"{workflow_geno}.gzi" if genome_zip else []
    log:
        f"{workflow_geno}.preprocess.log"
    params:
        genome_zip
    threads:
        workflow.cores
    conda:
        "envs/align.yaml"
    container:
        f"docker://pdimens/harpy:align_{VERSION}"
    shell: 
        """
        {{
            if (file {input} | grep -q compressed ) ;then
                # is regular gzipped, needs to be BGzipped
                seqtk seq {input} | bgzip -c > {output.geno}
            else
                cp -f {input} {output.geno}
            fi

            if [ "{params}" = "True" ]; then
                samtools faidx --gzi-idx {output.gzi} --fai-idx {output.fai} {output.geno}
            else
                samtools faidx --fai-idx {output.fai} {output.geno}
            fi

            minibwa index -t {threads} {output.geno} 
        }} 2> {log}
        """

rule align:
    input:
        multiext(workflow_geno, ".l2b", ".mbw"),
        ref   = workflow_geno,
        fastq = get_fq
    output:
        bam = temp("bwa/{sample}.bwa.bam"),
        tmp = temp(directory("bwa/{sample}_tmp"))
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
        mkdir -p {resources.tmpdir}
        {{
            minibwa map -t {threads} {params} {input.ref} {input.fastq} |
            samtools collate -T {resources.tmpdir} -O -u - 
        }} 2> {log} > {output.bam}
        """
