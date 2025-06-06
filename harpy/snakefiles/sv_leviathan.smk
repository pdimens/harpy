containerized: "docker://pdimens/harpy:latest"

import os
import re
import logging
from pathlib import Path

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake"]["log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+",
    population = r"[a-zA-Z0-9._-]+"

genomefile  = config["inputs"]["reference"]
bamlist     = config["inputs"]["alignments"]
bamdict     = dict(zip(bamlist,bamlist))
samplenames = {Path(i).stem for i in bamlist}
min_size      = config["min_size"]
min_bc      = config["min_barcodes"]
iterations  = config["iterations"]
small_thresh = config["variant_thresholds"]["small"]
medium_thresh = config["variant_thresholds"]["medium"]
large_thresh = config["variant_thresholds"]["large"]
duplcates_thresh = config["variant_thresholds"]["duplicates"]
extra       = config.get("extra", "") 
skip_reports = config["reports"]["skip"]
plot_contigs = config["reports"]["plot_contigs"]    
bn          = os.path.basename(genomefile)
if bn.lower().endswith(".gz"):
    workflow_geno = f"workflow/reference/{bn[:-3]}"
else:
    workflow_geno = f"workflow/reference/{bn}"

def get_alignments(wildcards):
    """returns a list with the bam file for the sample based on wildcards.sample"""
    r = re.compile(fr".*/({wildcards.sample})\.(bam|sam)$", flags = re.IGNORECASE)
    aln = list(filter(r.match, bamlist))
    return aln[0]

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
    shell:
        """
        samtools index {input} 2> {log}
        LRez index bam --threads {threads} -p -b {input} -o {output} 2>> {log}
        """

rule preprocess_reference:
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
    shell: 
        """
        seqtk seq {input} > {output.geno}
        samtools faidx --fai-idx {output.fai} {output.geno} 2> {log}
        bwa index {output.geno} 2>> {log}
        """

rule call_variants:
    input:
        bam = get_alignments,
        bc_idx = "lrez_index/{sample}.bci",
        genome = workflow_geno,
        genidx = multiext(workflow_geno, ".fai", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output:
        vcf = temp("vcf/{sample}.vcf"),
        candidates = "logs/leviathan/{sample}.candidates"
    log:  
        runlog = "logs/leviathan/{sample}.leviathan.log"
    params:
        min_size = f"-v {min_size}",
        min_bc = f"-c {min_bc}",
        iters  = f"-B {iterations}",
        small  = f"-s {small_thresh}",
        medium  = f"-m {medium_thresh}",
        large  = f"-l {large_thresh}",
        dupes  = f"-d {duplcates_thresh}",
        extra = extra
    threads:
        workflow.cores - 1
    conda:
        "envs/variants.yaml"
    shell:
        "LEVIATHAN -b {input.bam} -i {input.bc_idx} {params} -g {input.genome} -o {output.vcf} -t {threads} --candidates {output.candidates} 2> {log.runlog}"

rule sort_variants:
    priority: 100
    input:
        "vcf/{sample}.vcf"
    output:
        "vcf/{sample}.bcf"
    container:
        None
    shell:        
        "bcftools sort -Ob --output {output} {input} 2> /dev/null"

rule variant_stats:
    input: 
        "vcf/{sample}.bcf"
    output: 
        temp("reports/data/{sample}.sv.stats")
    container:
        None
    shell:
        """
        echo -e "sample\\tcontig\\tposition_start\\tposition_end\\tlength\\ttype\\tn_barcodes\\tn_pairs" > {output}
        bcftools query -f '{wildcards.sample}\\t%CHROM\\t%POS\\t%END\\t%SVLEN\\t%SVTYPE\\t%BARCODES\\t%PAIRS\\n' {input} >> {output}
        """

rule aggregate_variants:
    input:
        collect("reports/data/{sample}.sv.stats", sample = samplenames)
    output:
        "inversions.bedpe",
        "deletions.bedpe",
        "duplications.bedpe",
        "breakends.bedpe"
    run:
        with (
            open(output[0], "w") as inversions,
            open(output[1], "w") as deletions,
            open(output[2], "w") as duplications,
            open(output[3], "w") as breakends
        ):
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

rule configure_report:
    input:
        yaml = "workflow/report/_quarto.yml",
        scss = "workflow/report/_harpy.scss"
    output:
        yaml = temp("reports/_quarto.yml"),
        scss = temp("reports/_harpy.scss")
    run:
        import shutil
        for i,o in zip(input,output):
            shutil.copy(i,o)

rule sample_reports:
    input: 
        "reports/_quarto.yml",
        "reports/_harpy.scss",
        faidx     = f"{workflow_geno}.fai",
        statsfile = "reports/data/{sample}.sv.stats",
        qmd       = "workflow/report/leviathan.qmd"
    output:
        report = "reports/{sample}.leviathan.html",
        qmd = temp("reports/{sample}.leviathan.qmd")
    log:
        "logs/reports/{sample}.report.log"
    params:
        sample= lambda wc: "-P sample:" + wc.get('sample'),
        contigs= f"-P contigs:{plot_contigs}"
    conda:
        "envs/r.yaml"
    retries:
        3
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        FAIDX=$(realpath {input.faidx})
        STATS=$(realpath {input.statsfile})
        quarto render {output.qmd} --log {log} --quiet -P faidx:$FAIDX -P statsfile:$STATS {params}
        """

rule workflow_summary:
    default_target: True
    input: 
        vcf = collect("vcf/{sample}.bcf", sample = samplenames),
        bedpe_agg = collect("{sv}.bedpe", sv = ["inversions", "deletions","duplications", "breakends"]),
        reports = collect("reports/{sample}.leviathan.html", sample = samplenames) if not skip_reports else []
    params:
        min_size = f"-v {min_size}",
        min_bc = f"-c {min_bc}",
        iters  = f"-B {iterations}",
        extra = extra
    run:
        summary = ["The harpy sv leviathan workflow ran using these parameters:"]
        summary.append(f"The provided reference genome: {bn}")
        summary.append("The alignments were deconvolved using: leviathan_bx_shim.py")
        bc_idx = "The barcodes were indexed using:\n"
        bc_idx += "LRez index bam -p -b INPUT"
        summary.append(bc_idx)
        svcall = "Leviathan was called using:\n"
        svcall += f"\tLEVIATHAN -b INPUT -i INPUT.BCI -g GENOME {params}"
        summary.append(svcall)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['snakemake']['relative']}"
        summary.append(sm)
        with open("workflow/sv.leviathan.summary", "w") as f:
            f.write("\n\n".join(summary))