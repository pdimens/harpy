import os
from pathlib import Path
import re
from harpy.common.file_ops import pop_manifest

localrules: all, aggregate_variants
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+",

WORKFLOW   = config.get('Workflow', {})
PARAMETERS = config.get('Parameters', {})
INPUTS     = config['Inputs']
VERSION    = WORKFLOW.get('harpy-version', 'latest')

skip_reports  = WORKFLOW.get("reports", {}).get("skip", False)
plot_contigs  = WORKFLOW.get("reports", {}).get("plot-contigs", "default")
extra         = PARAMETERS.get("extra", "")
min_size      = PARAMETERS.get("min-size", 1000)
min_bc        = PARAMETERS.get("min-barcodes", 2)
iterations    = PARAMETERS.get("iterations", 50)
small_thresh  = PARAMETERS.get("variant-thresholds", {}).get("small", 95)
medium_thresh = PARAMETERS.get("variant-thresholds", {}).get("medium", 95)
large_thresh  = PARAMETERS.get("variant-thresholds", {}).get("large", 95)
dup_thresh    = PARAMETERS.get("variant-thresholds", {}).get("duplicates", 10)
genomefile    = INPUTS["reference"]
bamlist       = INPUTS["alignments"]
groupfile     = INPUTS.get("groupings", None)

plot_contigs  = ",".join(plot_contigs) if isinstance(plot_contigs, list) else plot_contigs
popdict       = pop_manifest(groupfile, bamlist) if groupfile else None
populations   = popdict.keys() if groupfile else None
target        = populations if groupfile else {Path(i).stem for i in bamlist}
bn            = os.path.basename(genomefile)
workflow_geno = f"workflow/reference/{bn[:-3]}" if bn.lower().endswith(".gz") else f"workflow/reference/{bn}"

def get_alignments(wildcards):
    """returns a list with the bam file for the sample based on wildcards.sample"""
    r = re.compile(fr".*/({wildcards.sample})\.(bam|sam)$", flags = re.IGNORECASE)
    aln = list(filter(r.match, bamlist))
    return aln[0]

rule process_reference:
    input:
        genomefile
    output: 
        multiext(workflow_geno, ".ann", ".bwt", ".pac", ".sa", ".amb"),
        geno = workflow_geno,
        fai  = f"{workflow_geno}.fai"
    log:
        f"{workflow_geno}.preprocess.log"
    conda:
        "envs/align.yaml"
    container:
        f"docker://pdimens/harpy:align_{VERSION}"
    shell: 
        """
        {{
            seqtk seq {input} > {output.geno}
            samtools faidx --fai-idx {output.fai} {output.geno}
            bwa index {output.geno}
        }} 2> {log}
        """

rule concat_groups:
    input: 
        bamfiles = lambda wc: collect("{samples}", samples = popdict[wc.sample]) 
    output:
        bam = temp("workflow/input/{sample}.bam"),
        bai = temp("workflow/input/{sample}.bam.bai")
    log:
        "logs/concat_groups/{sample}.concat.log"
    resources:
        mem_mb = 2000
    threads:
        workflow.cores
    shell:
        """
        {{
            djinn sam concat --bx {input} | 
            samtools sort -@ {threads} -O bam -l 0 -m {resources.mem_mb}M --write-index -o {output.bam}##idx##{output.bai}
        }} 2> {log}
        """

if popdict:
    rule index_barcodes:
        input: 
            bam = "workflow/input/{sample}.bam",
            bai = "workflow/input/{sample}.bam.bai"
        output:
            temp("lrez_index/{sample}.bci")
        log:
            "logs/lrez_index/{sample}.concat.log"
        threads:
            min(5, workflow.cores)
        conda:
            "envs/variants.yaml"
        container:
            f"docker://pdimens/harpy:variants_{VERSION}"
        shell:
            "LRez index bam -p -b {input.bam} -o {output} --threads {threads}"
else:
    rule index_barcodes:
        input: 
            get_alignments
        output:
            temp("lrez_index/{sample}.bci")
        log:
            "logs/process_alignments/{sample}.log"
        threads:
            min(10, workflow.cores)
        conda:
            "envs/variants.yaml"
        container:
            f"docker://pdimens/harpy:variants_{VERSION}"
        shell:
            """
            {{
                samtools index {input}
                LRez index bam --threads {threads} -p -b {input} -o {output}
            }} 2> {log}
            """

rule call_variants:
    input:
        bam    = "workflow/input/{sample}.bam" if popdict else get_alignments,
        bai    = "workflow/input/{sample}.bam.bai" if popdict else [],
        bc_idx = "lrez_index/{sample}.bci",
        genome = workflow_geno,
        genidx = multiext(workflow_geno, ".fai", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output:
        vcf = temp("vcf/{sample}.vcf"),
        candidates = "logs/leviathan/{sample}.candidates"
    log:  
        runlog = "logs/leviathan/{sample}.leviathan.log",
    params:
        min_size = f"-v {min_size}",
        min_bc = f"-c {min_bc}",
        iters  = f"-B {iterations}",
        small  = f"-s {small_thresh}",
        medium  = f"-m {medium_thresh}",
        large  = f"-l {large_thresh}",
        dupes  = f"-d {dup_thresh}",
        extra = extra
    threads:
        workflow.cores - 1
    conda:
        "envs/variants.yaml"
    container:
        f"docker://pdimens/harpy:variants_{VERSION}"
    shell:
        "LEVIATHAN -b {input.bam} -i {input.bc_idx} {params} -g {input.genome} -o {output.vcf} -t {threads} --candidates {output.candidates} 2> {log.runlog}"

rule sort_variants:
    priority: 100
    input:
        "vcf/{sample}.vcf"
    output:
        "vcf/{sample}.bcf"
    params:
        lambda wc: wc.sample
    shell:        
        "bcftools sort -Ob --output {output} {input} 2> /dev/null"

rule variant_stats:
    input: 
        "vcf/{sample}.bcf"
    output:
        temp("reports/data/{sample}.sv.stats")
    shell:
        """
        {{
            echo -e "sample\\tcontig\\tposition_start\\tposition_end\\tlength\\ttype\\tn_barcodes\\tn_pairs"
            bcftools query -f '{wildcards.sample}\\t%CHROM\\t%POS\\t%END\\t%SVLEN\\t%SVTYPE\\t%BARCODES\\t%PAIRS\\n' {input}
        }} > {output}
        """

rule aggregate_variants:
    input:
        collect("reports/data/{sample}.sv.stats", sample = target)
    output:
        "inversions.bedpe",
        "deletions.bedpe",
        "duplications.bedpe",
        "breakends.bedpe"
    run:
        with open(output[0], "w") as inversions, open(output[1], "w") as deletions, open(output[2], "w") as duplications, open(output[3], "w") as breakends:
            header = ["sample","contig","position_start","position_end","length","type","n_barcodes","n_pairs"]
            _ = inversions.write("\t".join(header) + "\n")
            _ = deletions.write("\t".join(header) + "\n")
            _ = duplications.write("\t".join(header) + "\n")
            _ = breakends.write("\t".join(header) + "\n")
            for varfile in input:
                with open(varfile, "r") as f_in:
                    # skip header
                    f_in.readline()
                    while True:
                        line = f_in.readline()
                        if not line:
                            break
                        record = line.rstrip().split("\t")
                        if record[5] == "INV":
                            _ = inversions.write(line)
                        elif record[5] == "DEL":
                            _ = deletions.write(line)
                        elif record[5] == "DUP":
                            _ = duplications.write(line)
                        elif record[5] == "BND":
                            _ = breakends.write(line)

rule report:
    input: 
        faidx = f"{workflow_geno}.fai",
        stats = collect("{var}.bedpe", var = ['inversions', 'deletions', 'duplications', 'breakends']),
        ipynb = "workflow/sv.ipynb"
    output:
        tmp = temp("reports/leviathan.summary.tmp.ipynb"),
        ipynb = "reports/leviathan.summary.ipynb"
    log:
        "logs/report.log"
    params:
        f"-p indir {os.getcwd()}",
        f"-p faidx " + os.path.abspath(f"{workflow_geno}.fai"),
        f"-p contigs {plot_contigs}" if plot_contigs != "default" else ""
    shell:
        """
        {{
            papermill -k xpython --no-progress-bar --log-level ERROR {input.ipynb} {output.tmp} {params}
            harpy-utils process-notebook {output.tmp} LEVIATHAN
        }} 2> {log} > {output.ipynb}
        """

rule all:
    default_target: True
    input:
        vcf = collect("vcf/{sample}.bcf", sample = target),
        bedpe_agg = collect("{sv}.bedpe", sv = ["inversions", "deletions","duplications", "breakends"]),
        agg_report = "reports/leviathan.summary.ipynb" if not skip_reports else []
