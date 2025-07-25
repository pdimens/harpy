containerized: "docker://pdimens/harpy:latest"

import os
import logging
from pathlib import Path

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake"]["log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

ploidy 		= config["ploidy"]
extra 	    = config.get("extra", "") 
regions_input = config["inputs"]["regions"]
skip_reports = config["reports"]["skip"]
bamlist     = config["inputs"]["alignments"]
bamdict     = dict(zip(bamlist, bamlist))
genomefile 	= config["inputs"]["reference"]
bn          = os.path.basename(genomefile)
if bn.lower().endswith(".gz"):
    genome_zip  = True
    bn = bn[:-3]
else:
    genome_zip  = False
workflow_geno = f"workflow/reference/{bn}"
groupings 	= config["inputs"].get("groupings", [])
samplenames = {Path(i).stem for i in bamlist}
sampldict = dict(zip(bamlist, samplenames))

if os.path.exists(regions_input):
    with open(regions_input, "r") as reg_in:
        intervals = set()
        for line in reg_in:
            cont,startpos,endpos = line.split()
            intervals.add(f"{cont}:{max(int(startpos),1)}-{int(endpos)}")
    regions = dict(zip(intervals, intervals))
else:
    intervals = [regions_input]
    regions = {f"{regions_input}" : f"{regions_input}"}

rule preprocess_groups:
    input:
        grp = groupings
    output:
        grp = "workflow/sample.groups"
    run:
        with open(input.grp, "r") as infile, open(output.grp, "w") as outfile:
            _ = [outfile.write(i) for i in infile.readlines() if not i.lstrip().startswith("#")]

rule preprocess_reference:
    input:
        genomefile
    output: 
        geno = workflow_geno,
        fai = f"{workflow_geno}.fai"
    log:
        f"{workflow_geno}.preprocess.log"
    container:
        None
    shell: 
        """
        seqtk seq {input} > {output}
        samtools faidx --fai-idx {output.fai} {output.geno} 2> {log}
        """

rule index_alignments:
    input:
        lambda wc: bamdict[wc.bam]
    output:
        "{bam}.bai"
    container:
        None
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
    shell:
        """
        freebayes -f {input.reference} -L {input.bamlist} {params} 2> {log} |
            bcftools sort - --output {output.bcf} --write-index 2> /dev/null
        """

rule concat_list:
    input:
        bcfs = collect("regions/{part}.bcf", part = intervals),
    output:
        bcf = "logs/bcf.files"
    run:
        with open(output.bcf, "w") as fout:
            for bcf in input.bcfs:
                _ = fout.write(f"{bcf}\n")

rule concat_variants:
    input:
        collect("regions/{part}.{ext}", part = intervals, ext = ["bcf", "bcf.csi"]),
        filelist = "logs/bcf.files"
    output:
        temp("variants.raw.unsort.bcf")
    log:
        "logs/concat.log"
    threads:
        workflow.cores
    container:
        None
    shell:  
        "bcftools concat -f {input.filelist} --threads {threads} --naive -Ob -o {output} 2> {log}"

rule sort_variants:
    input:
        "variants.raw.unsort.bcf"
    output:
        bcf = "variants.raw.bcf",
        csi = "variants.raw.bcf.csi"
    container:
        None
    shell:
        "bcftools sort --write-index -Ob -o {output.bcf} {input} 2> /dev/null"

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
    threads:
        workflow.cores
    container:
        None
    shell:
        "bcftools norm --threads {threads} {params} -o {output.bcf} -f {input.genome} {input.bcf} 2> {log}"    

rule general_stats:
    input:
        genome  = workflow_geno,
        ref_idx = f"{workflow_geno}.fai",
        bcf     = "variants.{type}.bcf",
        idx     = "variants.{type}.bcf.csi"
    output:
        "reports/data/variants.{type}.stats",
    container:
        None
    shell:
        """
        bcftools stats -s "-" --fasta-ref {input.genome} {input.bcf} > {output} 2> /dev/null
        """

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

rule variant_report:
    input: 
        "reports/_quarto.yml",
        "reports/_harpy.scss",
        data = "reports/data/variants.{type}.stats",
        qmd  = "workflow/report/bcftools_stats.qmd"
    output:
        report = "reports/variants.{type}.html",
        qmd = temp("reports/variants.{type}.qmd")
    params:
        lambda wc: "-P vcf:variants." + wc.get("type")
    log:
        "logs/variants.{type}.report.log"
    conda:
        "envs/r.yaml"
    retries:
        3
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        INPATH=$(realpath {input.data})
        quarto render {output.qmd} --log {log} --quiet -P infile:$INPATH {params}
        """

rule workflow_summary:
    default_target: True
    input:
        vcf = collect("variants.{file}.bcf", file = ["raw","normalized"]),
        reports = collect("reports/variants.{file}.html", file = ["raw","normalized"]) if not skip_reports else []
    params:
        ploidy = f"-p {ploidy}",
        populations = f"--populations {groupings}" if groupings else '',
        extra = extra
    run:
        summary = ["The harpy snp freebayes workflow ran using these parameters:"]
        summary.append(f"The provided reference genome: {bn}")
        summary.append(f"Genomic positions for which variants were called: {regions_input}")
        varcall = "The freebayes parameters:\n"
        varcall += f"\tfreebayes -f REFERENCE -L samples.list -r REGION {params} |\n"
        varcall += f"\tbcftools sort -"
        summary.append(varcall)
        merged = "The variants identified in the intervals were merged into the final variant file using:\n"
        merged += "\tbcftools concat -f bcf.files -a --remove-duplicates"
        summary.append(merged)
        normalize = "The variants were normalized using:\n"
        normalize += "\tbcftools norm -m -both -d both -c w"
        summary.append(normalize)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['snakemake']['relative']}"
        summary.append(sm)
        with open("workflow/snp.freebayes.summary", "w") as f:
            f.write("\n\n".join(summary))
