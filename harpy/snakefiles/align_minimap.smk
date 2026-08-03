import os
import re

wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

WORKFLOW    = config.get('Workflow') or {}
PARAMETERS = config.get('Parameters') or {}
REPORTS    = WORKFLOW.get("reports") or {} 
INPUTS     = config['Inputs']
VERSION    = WORKFLOW.get('harpy-version', 'latest')

tech         = PARAMETERS.get("technology", "sr") 
illumina_old = PARAMETERS.get("illumina-format-old", False)
extra 		 = PARAMETERS.get("extra", "") 
fqlist       = INPUTS["fastq"]
genomefile 	 = INPUTS["reference"]

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
        idx = multiext(workflow_geno, ".mmi"),
        fai = f"{workflow_geno}.fai",
        gzi = f"{workflow_geno}.gzi" if genome_zip else []
    log:
        f"{workflow_geno}.preprocess.log"
    params:
        genome_zip
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
            minimap2 -d {output.idx} {output.geno} 
        }} 2> {log}
        """

rule align:
    input:
        ref = workflow_geno + ".mmi",
        fastq = get_fq
    output:
        bam = temp("minimap/{sample}.minimap.bam"),
        tmp = temp(directory("minimap/{sample}_tmp"))
    log:
        "logs/minimap/{sample}.minimap.log"
    params:
        RG_tag = lambda wc: "-R \"@RG\\tID:" + wc.get("sample") + "\\tSM:" + wc.get("sample") + "\"",
        tech = f"-ax map-{technology}",
        static = "--MD -y" if illumina_old else "--MD",
        extra = extra
    threads:
        12
    resources:
        tmpdir = lambda wc: f"minimap/{wc.sample}_tmp"
    conda:
        "envs/align.yaml"
    container:
        f"docker://pdimens/harpy:align_{VERSION}"
    shell:
        """
        mkdir -p {resources.tmpdir}
        {{
            minimap2 -t {threads} {params} {input} |
            samtools collate -T {resources.tmpdir} -O -u - 
        }} 2> {log} > {output.bam}
        """
