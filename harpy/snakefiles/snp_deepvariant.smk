import os
from pathlib import Path

localrules: all
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

WORKFLOW   = config.get('Workflow') or {}
PARAMETERS = config.get('Parameters') or {}
REPORTS    = WORKFLOW.get("reports") or {} 
INPUTS     = config['Inputs']
VERSION    = WORKFLOW.get('harpy-version', 'latest')

skip_reports  = REPORTS.get("skip", False)
ploidy 		  = PARAMETERS.get("ploidy", 2)
extra 	      = PARAMETERS.get("extra", "") 
bamlist       = INPUTS["alignments"]
genomefile 	  = INPUTS["reference"]
region_input = INPUTS["regions"]
keep_invar   = PARAMETERS.get("keep-invariant", False)
gpu          = WORKFLOW.get('gpu', False)

# attempt to get processed, then source, then nothing
grp          = INPUTS.get("groupings") or {}
if grp:
    groupings = grp.get("processed", [])
    if isinstance(groupings, str) and not os.path.isfile(groupings):
        groupings = grp.get("source", [])
else:
    groupings = []

bamdict       = dict(zip(bamlist, bamlist))
samplenames   = {Path(i).stem for i in bamlist}
sampldict     = dict(zip(bamlist, samplenames))
bn                = os.path.basename(genomefile)
genome_zip        = True if bn.lower().endswith(".gz") else False
workflow_geno     = f"workflow/reference/{bn}"
workflow_geno_idx = f"{workflow_geno}.gzi" if genome_zip else f"{workflow_geno}.fai"

if os.path.exists(region_input):
    with open(region_input, "r") as reg_in:
        intervals = set()
        for line in reg_in:
            cont,startpos,endpos = line.split()
            intervals.add(f"{cont}:{max(int(startpos),1)}-{int(endpos)}")
    regions = dict(zip(intervals, intervals))
else:
    intervals = [region_input]
    regions   = {f"{region_input}" : f"{region_input}"}

def get_alignments(wildcards):
    """returns a list with the bam file for the sample based on wildcards.sample"""
    r = re.compile(fr".*/({wildcards.sample})\.(bam|sam)$", flags = re.IGNORECASE)
    aln = list(filter(r.match, bamlist))
    return aln[0]

def get_align_index(wildcards):
    """returns a list with the bai index file for the sample based on wildcards.sample"""
    r = re.compile(fr"(.*/{wildcards.sample})\.(bam|sam)$", flags = re.IGNORECASE)
    aln = list(filter(r.match, bamlist))
    return aln[0] + ".bai"

rule process_reference:
    input:
        genomefile
    output: 
        geno = workflow_geno,
        fai = f"{workflow_geno}.fai",
        gzi = f"{workflow_geno}.gzi" if genome_zip else []
    log:
        f"{workflow_geno}.preprocess.log"
    params:
        f"--gzi-idx {workflow_geno}.gzi" if genome_zip else ""
    shell: 
        """
        {{
            if (file {input} | grep -q compressed ) ;then
                # is regular gzipped, needs to be BGzipped
                seqtk seq {input} | bgzip -c > {output.geno}
            else
                ln -s {input} {output.geno}
            fi
            samtools faidx {params} --fai-idx {output.fai} {output.geno}
        }} 2> {log}
        """

rule index_alignments:
    input:
        lambda wc: bamdict[wc.bam]
    output:
        "{bam}.bai"
    shell:
        "samtools index {input}"

# either vcf or gvcf
rule call_variants:
    input:
        get_alignments_index,
        bam = get_alignments,
        f"{workflow_geno}.fai",
        reference = workflow_geno
    output:
        dir("deepvariant/{sample}"),
        vcf = temp("samples/{sample}.vcf")
    log:
        "logs/{sample}.deepvariant.log"
    params:
        "--model_type=WGS",
        "--use_gpu" if gpu else "",
        lambda wc: "--intermediate_results_dir=deepvariant/{wc.sample}",
        extra = extra,
        "--output_gvcf=" if keep_invar else "--output_vcf="
    threads:
        4
    container:
        "docker://google/deepvariant:1.10.0"
    shell:
        """
        mkdir -p deepvariant/{wildcards.sample}
        run_deepvariant --ref={input.reference} --reads={input.bam} --num_shards={threads} {params}{output.vcf} &> {log}
        """

rule sort_variants:
    input:
        "samples/{sample}.vcf"
    output:
        temp("samples/{sample}.bcf.csi"),
        bcf = temp("samples/{sample}.bcf")
    threads:
        2
    shell:
        "bcftools sort -Ou --write-index {input} > {output.bcf}"


rule concat_samples:
    input:
        collect("samples/{sample}.bcf.csi", part = samplenames),
        bcf = collect("samples/{sample}.bcf", part = samplenames)
    output:
        concatlist = temp("logs/bcf.files"),
        bcf = "variants.raw.bcf",
        csi = "variants.raw.bcf.csi"
    log:
        "logs/concat_sort.log"
    threads:
        workflow.cores
    params:
        workflow.cores - 1
    resources:
        mem_mb = 8000
    shell:  
        """
        printf '%s\\n' {input.bcf} > {output.concatlist}
        {{
            bcftools merge -@ {params} --no-version -f {output.concatlist} --max-mem {resources}M |
            bcftools sort - --write-index -Ob -o {output.bcf}
        }} 2> {log}
        """

rule realign_indels:
    input:
        genome  = workflow_geno,
        bcf     = "variants.raw.bcf",
        idx     = "variants.raw.bcf.csi"
    output:
        bcf = "variants.normalized.bcf",
        idx = "variants.normalized.bcf.csi"
    log:
        "logs/variants.normalized.log"
    threads:
        workflow.cores
    params:
        "-m -both -d both --write-index -Ob -c w"
    shell:
        "bcftools norm --threads {threads} {params} -o {output.bcf} -f {input.genome} {input.bcf} 2> {log}"    

rule variant_report:
    input: 
        genome  = workflow_geno,
        ref_idx = f"{workflow_geno}.fai",
        bcf     = "variants.{type}.bcf",
        idx     = "variants.{type}.bcf.csi",
        ipynb  = "workflow/bcftools_stats.ipynb"
    output:
        data = temp("reports/data/variants.{type}.stats"),
        tmp =  temp("reports/variants.{type}.tmp.ipynb"),
        ipynb = "reports/variants.{type}.ipynb"
    log:
        "logs/variants.{type}.report.log"
    params:
        lambda wc: "-p infile " + os.path.abspath(f"reports/data/variants.{wc.type}.stats")
    shell:
        """
        export IPYTHONDIR=/tmp/ipython-snp-{wildcards.type}
        {{
            bcftools stats -s "-" --fasta-ref {input.genome} {input.bcf} > {output.data} 
            papermill -k ipython-harpy --no-progress-bar --log-level ERROR {input.ipynb} {output.tmp} {params}
            harpy-utils process-notebook {output.tmp} "Variants ({wildcards.type})" > {output.ipynb}
        }} 2> {log}
        """

rule all:
    default_target: True
    input:
        vcf = collect("variants.{file}.bcf", file = ["raw","normalized"]),
        reports = collect("reports/variants.{file}.ipynb", file = ["raw","normalized"]) if not skip_reports else []

