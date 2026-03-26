import os
from pathlib import Path

localrules: all
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

WORKFLOW   = config.get('Workflow', {})
PARAMETERS = config.get('Parameters', {})
INPUTS     = config['Inputs']
VERSION    = WORKFLOW.get('harpy-version', 'latest')

skip_reports  = WORKFLOW.get("reports", {}).get("skip", False)
ploidy 		  = PARAMETERS.get("ploidy", 2)
extra 	      = PARAMETERS.get("extra", "") 
bamlist       = INPUTS["alignments"]
genomefile 	  = INPUTS["reference"]
regions_input = INPUTS["regions"]
groupings 	  = INPUTS.get("groupings", [])

bamdict       = dict(zip(bamlist, bamlist))
bn            = os.path.basename(genomefile)
genome_zip    = bn.lower().endswith(".gz")
bn            = bn[:-3] if genome_zip else bn
workflow_geno = f"workflow/reference/{bn}"
samplenames   = {Path(i).stem for i in bamlist}
sampldict     = dict(zip(bamlist, samplenames))

if os.path.exists(regions_input):
    with open(regions_input, "r") as reg_in:
        intervals = set()
        for line in reg_in:
            cont,startpos,endpos = line.split()
            intervals.add(f"{cont}:{max(int(startpos),1)}-{int(endpos)}")
    regions = dict(zip(intervals, intervals))
else:
    intervals = [regions_input]
    regions   = {f"{regions_input}" : f"{regions_input}"}

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

rule index_alignments:
    input:
        lambda wc: bamdict[wc.bam]
    output:
        "{bam}.bai"
    shell:
        "samtools index {input}"

rule bam_list:
    input: 
        collect("{bam}.bai", bam = bamlist),
        bam = bamlist
    output:
        smp = "workflow/freebayes.input"
    run:
        with open(output.smp, "w") as fout:
            for bamfile in input.bam:
                _ = fout.write(bamfile + "\n")

rule call_variants:
    input:
        bamlist,
        collect("{bam}.bai", bam = bamlist),
        "workflow/sample.groups" if groupings else [],
        f"{workflow_geno}.fai",
        reference = workflow_geno,
        bamlist  = "workflow/freebayes.input"
    output:
        bcf = temp("regions/{part}.bcf"),
        idx = temp("regions/{part}.bcf.csi")
    log:
        "logs/{part}.freebayes.log"
    params:
        region = lambda wc: "-r " + regions[wc.part],
        ploidy = f"-p {ploidy}",
        static = "--strict-vcf",
        populations = "--populations workflow/sample.groups" if groupings else "",
        extra = extra
    conda:
        "envs/variants.yaml"
    container:
        f"docker://pdimens/harpy:variants_{VERSION}"
    shell:
        """
        freebayes -f {input.reference} -L {input.bamlist} {params} 2> {log} |
            bcftools sort - --output {output.bcf} --write-index 2> /dev/null
        """

rule concat_variants:
    input:
        collect("regions/{part}.bcf.csi", part = intervals),
        bcf = collect("regions/{part}.bcf", part = intervals)
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
    shell:  
        """
        for i in {input.bcf}; do echo $i; done >> {output.concatlist}
        {{
            bcftools concat -f {output.concatlist} --threads {params} --naive |
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
        {{
            bcftools stats -s "-" --fasta-ref {input.genome} {input.bcf} > {output.data} 
            papermill -k xpython --no-progress-bar --log-level ERROR {input.ipynb} {output.tmp} {params}
            harpy-utils process-notebook {output.tmp} variants.{wildcards.type}
        }} 2> {log} > {output.ipynb}
        """

rule all:
    default_target: True
    input:
        vcf = collect("variants.{file}.bcf", file = ["raw","normalized"]),
        reports = collect("reports/variants.{file}.ipynb", file = ["raw","normalized"]) if not skip_reports else []

