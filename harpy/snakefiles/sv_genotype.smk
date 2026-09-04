
localrules: all
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

WORKFLOW   = config.get('Workflow') or {}
PARAMETERS = config.get('Parameters') or {}
REPORTS    = WORKFLOW.get("reports") or {} 
INPUTS     = config['Inputs']
VERSION    = WORKFLOW.get('harpy-version', 'latest')

reference  = INPUTS["reference"]
fqlist     = INPUTS["fastq"]
vcf        = INPUTS["vcf"]

regionsize = PARAMETERS.get("region-size", 10000)
inaccuracy = PARAMETERS.get("inaccuracy", 0)
mindiff    = PARAMETERS.get("likelihood-min-diff", 20)
proberror  = PARAMETERS.get("likelihood-prob-error", [0.2,0.1,0.02,0.008])

bn_r          = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
samplenames   = {re.sub(bn_r, "", os.path.basename(i), flags = re.IGNORECASE) for i in fqlist}
d             = dict(zip(samplenames, samplenames))

def get_fq(wildcards):
    # returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    r = re.compile(fr".*/({re.escape(wildcards.sample)}){bn_r}", flags = re.IGNORECASE)
    return sorted(list(filter(r.match, fqlist))[:2])

rule svjeditag:
    input:
        fastq = get_fq,
        reference = reference,
        vcf = vcf
    output:
        "genotypes.vcf"
    params:
        f"--regionSize {regionsize}",
        f"--inaccuracy {inaccuracy}",           
        f"--likelihood_MinDiff {mindiff}",   
        f"--likelihood_ProbError {proberror}",
        lambda wc: f"-o {wc.sample}.genotypes"
    threads:
        workflow.cores
    conda:
        "envs/variants.yaml"
    container:
        f"docker://pdimens/harpy:variants_{VERSION}"
    shell:
        "svjedi-tag -v {input.vcf} -r {input.reference} -q {input.fastq} {params}"

rule all:
    default_target: True
    input:
        collect("{sample}.genotypes.vcf", sample = samplenames),
