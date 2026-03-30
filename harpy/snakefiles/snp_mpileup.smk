import os
from pathlib import Path

localrules: all, concat_logs
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

WORKFLOW   = config.get('Workflow') or {}
PARAMETERS = config.get('Parameters') or {}
REPORTS    = WORKFLOW.get("reports") or {} 
INPUTS     = config['Inputs']
VERSION    = WORKFLOW.get('harpy-version', 'latest')

skip_reports = REPORTS.get("skip", False)
ploidy 		 = PARAMETERS.get("ploidy", 2)
mp_extra 	 = PARAMETERS.get("extra", "")
bamlist      = INPUTS["alignments"]
genomefile 	 = INPUTS["reference"]
region_input = INPUTS["regions"]
# attempt to get processed, then source, then nothing
grp          = INPUTS.get("groupings") or {}
if grp:
    groupings = grp.get("processed", [])
    if not os.path.isfile(groupings):
        groupings.get("source") or []
else:
    groupings = []

bamdict           = dict(zip(bamlist, bamlist))
bn                = os.path.basename(genomefile)
genome_zip        = True if bn.lower().endswith(".gz") else False
workflow_geno     = f"workflow/reference/{bn}"
workflow_geno_idx = f"{workflow_geno}.gzi" if genome_zip else f"{workflow_geno}.fai"
samplenames       = {Path(i).stem for i in bamlist}

if os.path.isfile(region_input):
    with open(region_input, "r") as reg_in:
        intervals = set()
        for line in reg_in:
            cont,startpos,endpos = line.split()
            intervals.add(f"{cont}:{max(int(startpos),1)}-{int(endpos)}")
    regions = dict(zip(intervals, intervals))
else:
    intervals = [region_input]
    regions   = {f"{region_input}" : f"{region_input}"}

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
                cp -f {input} {output.geno}
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

rule bam_list:
    input: 
        collect("{bam}.bai", bam = bamlist),
        bam = bamlist
    output:
        temp("workflow/mpileup.input")
    run:
        with open(output[0], "w") as fout:
            for bamfile in input.bam:
                _ = fout.write(bamfile + "\n")

rule call_genotypes:
    input:
        bamlist,
        collect("{bam}.bai", bam = bamlist),
        f"{workflow_geno}.fai",
        "workflow/sample.groups" if groupings else [],
        bamlist = "workflow/mpileup.input",
        genome  = workflow_geno,
    output: 
        vcf = temp("call/{part}.vcf"),
        logfile = temp("logs/mpileup/{part}.mpileup.log")
    log:
        "logs/call/{part}.call.log"
    params:
        region = lambda wc: "-r " + regions[wc.part],
        annot_mp = "-a AD,INFO/FS",
        extra = mp_extra,
        ploidy = f"--ploidy {ploidy}",
        annot_call = "-a GQ,GP",
        groups = "--group-samples workflow/sample.groups" if groupings else "--group-samples -"
    shell:
        """
        bcftools mpileup --threads {threads} --fasta-ref {input.genome} --bam-list {input.bamlist} -Ou {params.region} {params.annot_mp} {params.extra} 2> {output.logfile} |
            bcftools call -o {output.vcf} --multiallelic-caller --variants-only {params.ploidy} {params.annot_call} {params.groups} 2> {log}
        """

rule sort_variants:
    input:
        bcf = temp("call/{part}.vcf")
    output:
        bcf = temp("sort/{part}.bcf"),
        idx = temp("sort/{part}.bcf.csi")
    log:
        "logs/sort/{part}.sort.log"
    shell:
        "bcftools sort --output {output.bcf} --write-index {input.bcf} 2> {log}"

rule concat_variants:
    input:
        collect("sort/{part}.bcf.csi", part = intervals),
        bcf = collect("sort/{part}.bcf", part = intervals)
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

rule concat_logs:
    input:
        collect("logs/mpileup/{part}.mpileup.log", part = intervals)
    output:
        "logs/mpileup.log"
    shell:  
        """
        for i in {input}; do
            interval=$(basename "$i" .mpileup.log)
            awk -v prefix="$interval" '{{print prefix "\t" $0}}' "$i"
        done >> {output}
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
    params:
        "-m -both -d both --write-index -Ob -c w"
    threads:
        workflow.cores
    shell:
        "bcftools norm --threads {threads} {params} -o {output.bcf} -f {input.genome} {input.bcf} 2> {log}"

rule variant_report:
    input: 
        genome  = workflow_geno,
        bcf     = "variants.{type}.bcf",
        idx     = "variants.{type}.bcf.csi",
        ipynb  = "workflow/bcftools_stats.ipynb"
    output:
        data = temp("reports/data/variants.{type}.stats"),
        tmp = temp("reports/variants.{type}.tmp.ipynb"),
        ipynb = temp("reports/variants.{type}.ipynb")
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
        vcf = collect("variants.{file}.bcf", file = ["raw", "normalized"]),
        agg_log = "logs/mpileup.log",
        reports = collect("reports/variants.{file}.ipynb", file = ["raw", "normalized"]) if not skip_reports else []
