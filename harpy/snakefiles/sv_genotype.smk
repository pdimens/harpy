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

def get_fq1(wildcards):
    '''returns a list of fastq files for read 1 based on *wildcards.sample*'''
    r = re.compile(fr"(.*/{re.escape(wildcards.sample)})([_\.]1|[_\.]F|[_\.]R1(?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$", flags = re.IGNORECASE)
    return list(filter(r.match, fqlist))

def get_fq2(wildcards):
    '''returns a list of fastq files for read 2 based on *wildcards.sample*'''
    r = re.compile(fr"(.*/{re.escape(wildcards.sample)})([_\.]2|[_\.]R|[_\.]R2(?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$", flags = re.IGNORECASE)
    return list(filter(r.match, fqlist))

rule construct_graph:
    input:
        vcf = vcf,
        ref = reference
    output:
        "graph/graph.gfa"
    log:
        "logs/graph.build.log"
    conda:
        "envs/variants.yaml"
    container:
        f"docker://pdimens/harpy:variants_{VERSION}"
    shell:
        #TODO SCRIPT DIR, FIGURE THAT OUT
        "python3 script_dir/construct_graph.py -v {input.vcf} -r {input.reference} -o {output} 2> {log}"

rule index_graph:
    input:
        "graph/graph.gfa"
    output:
        "graph/graph.gaf",
        "graph/graph.giraffe.gbz",
        "graph/graph.shortread.withzip.min",
        "graph/graph.shortread.zipcodes",
        "graph/graph.dist"
    log:
        "logs/graph.index.log"
    conda:
        "envs/variants.yaml"
    container:
        f"docker://pdimens/harpy:variants_{VERSION}"
    shell:
        "vg autoindex --workflow sr-giraffe -g {input} -p graph/graph 2> {log}"

rule map_to_graph:
    input:
        R1 = get_fq1,
        R2 = get_fq2,
        gbz  = "graph/graph.giraffe.gbz",
        zmin  = "graph/graph.shortread.withzip.min",
        zipc  = "graph/graph.shortread.zipcodes",
        dist = "graph/graph.dist",
    output:
        "map/{sample}.gaf"
    threads:
        6
    conda:
        "envs/variants.yaml"
    container:
        f"docker://pdimens/harpy:variants_{VERSION}"
    shell:
        "vg giraffe -t {threads} -Z {input.gbz} -m {input.zmin} -z {input.zipc} -d {input.dist} -f {input.R1} -f {input.R2} "
        "-o gaf --named-coordinates --comments-as-tags > {output} 2> {log}"

rule call_genotypes:
    input:
        aln = "map/{sample}.gaf",
        vcf = vcf,
        gfa = "graph/graph.gfa"
    output:
        "{sample}.genotypes.vcf"
    params:
        f"-s {regionsize}",
        f"-i {inaccuracy}",
        f"-d {mindiff}",
        f"-e {proberror[0]} {proberror[1]} {proberror[2]} {proberror[3]}"
    shell:
        "python3 scriptdir/predict_genotype.py -a {input.aln} -v {input.vcf} -o {output} -g {input.gfa} {params}"

rule all:
    default_target: True
    input:
        collect("{sample}.genotypes.vcf", sample = samplenames)
