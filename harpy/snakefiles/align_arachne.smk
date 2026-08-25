import os
import re

localrules: all
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

WORKFLOW    = config.get('Workflow') or {}
PARAMETERS = config.get('Parameters') or {}
REPORTS    = WORKFLOW.get("reports") or {} 
INPUTS     = config['Inputs']
VERSION    = WORKFLOW.get('harpy-version', 'latest')

molecule_distance = PARAMETERS.get("distance-threshold", 0)
fqlist            = INPUTS["fastq"]
genomefile 	      = INPUTS["reference"]
centromeres 	  = INPUTS.get("centromeres", None)

# overwrite what's in align.smk
bx_tag            = True 
vx_tag            = True
ignore_bx         = False

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
        multiext(workflow_geno, '.amb', '.ann', '.bwt', '.pac', '.sa'),
        geno = workflow_geno,
        fai = f"{workflow_geno}.fai",
        gzi = f"{workflow_geno}.gzi" if genome_zip else []
    log:
        f"{workflow_geno}.preprocess.log"
    params:
        genome_zip
    #conda:
    #    "envs/align.yaml"
    #container:
    #    f"docker://pdimens/harpy:align_{VERSION}"
    shell: 
        """
        {{
            if (file {input} | grep -q compressed ) ;then
                # is regular gzipped, needs to be BGzipped
                seqtk seq {input} | bgzip -c > {output.geno}
            else
                ln -s {input} {output.geno}
            fi

            if [ "{params}" = "True" ]; then
                samtools faidx --gzi-idx {output.gzi} --fai-idx {output.fai} {output.geno}
            else
                samtools faidx --fai-idx {output.fai} {output.geno}
            fi

            arachne index {output.geno}
        }} 2> {log}
        """

rule arachne_prep:
    input:
        get_fq
    output:
        R1 = temp("arachne-prep/{sample}.R1.fq.gz"),
        R2 = temp("arachne-prep/{sample}.R2.fq.gz"),
    threads:
        4
    log:
        "logs/preprocess/{sample}.prep.log"
#    conda:
#        "envs/align.yaml"
#    container:
#        f"docker://pdimens/harpy:align_{VERSION}"
    shell:
        """
        arachne prep -t {threads} arachne-prep/{wildcards.sample} {input} 2> {log}
        mv arachne-prep/{wildcards.sample}.arachne.R1.fq.gz {output.R1}
        mv arachne-prep/{wildcards.sample}.arachne.R2.fq.gz {output.R2}
        """

rule align:
    input:
        multiext(workflow_geno, '.amb', '.ann', '.bwt', '.pac', '.sa'),
        ref   = workflow_geno,
        R1 = "arachne-prep/{sample}.R1.fq.gz",
        R2 = "arachne-prep/{sample}.R2.fq.gz",
        centromeres = centromeres if centromeres else []
    output:
        bam = temp("arachne/{sample}.arachne.bam"),
        tmp = temp(directory("arachne/{sample}_tmp"))
    log:
        "logs/arachne/{sample}.arachne.log"
    params:
        f"-d {molecule_distance}",
        lambda wc: "-s " + wc.get("sample")
    threads:
        12
    resources:
        tmpdir = lambda wc: f"arachne/{wc.sample}_tmp"
#    conda:
#        "envs/align.yaml"
#    container:
#        f"docker://pdimens/harpy:align_{VERSION}"
    shell:
        """
        {{
            mkdir -p {resources.tmpdir}
            arachne align -t {threads} {params} {input.ref} {input.R1} {input.R2} |
            samtools collate -T {resources.tmpdir} -O -u - 
        }} > {output.bam} 2> {log}
        """
