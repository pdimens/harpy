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
mp_extra 	= config.get("extra", "")
skip_reports = config["reports"]["skip"]
bamlist     = config["inputs"]["alignments"]
bamdict     = dict(zip(bamlist, bamlist))
genomefile 	= config["inputs"]["reference"]
bn          = os.path.basename(genomefile)
workflow_geno = f"workflow/reference/{bn}"
genome_zip  = True if bn.lower().endswith(".gz") else False
workflow_geno_idx = f"{workflow_geno}.gzi" if genome_zip else f"{workflow_geno}.fai"
groupings 	= config["inputs"].get("groupings", [])
region_input = config["inputs"]["regions"]
samplenames = {Path(i).stem for i in bamlist}

if os.path.isfile(region_input):
    with open(region_input, "r") as reg_in:
        intervals = set()
        for line in reg_in:
            cont,startpos,endpos = line.split()
            intervals.add(f"{cont}:{max(int(startpos),1)}-{int(endpos)}")
    regions = dict(zip(intervals, intervals))
else:
    intervals = [region_input]
    regions = {f"{region_input}" : f"{region_input}"}

rule preprocess_reference:
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
    container:
        None
    shell: 
        """
        if (file {input} | grep -q compressed ) ;then
            # is regular gzipped, needs to be BGzipped
            seqtk seq {input} | bgzip -c > {output.geno}
        else
            cp -f {input} {output.geno}
        fi
        samtools faidx {params} --fai-idx {output.fai} {output.geno} 2> {log}
        """

rule preprocess_groups:
    input:
        grp = groupings
    output:
        grp = "workflow/sample.groups"
    run:
        with open(input.grp, "r") as infile, open(output.grp, "w") as outfile:
            _ = [outfile.write(i) for i in infile.readlines() if not i.lstrip().startswith("#")]

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
    threads:
        1
    container:
        None
    shell:
        """
        bcftools mpileup --threads {threads} --fasta-ref {input.genome} --bam-list {input.bamlist} -Ou {params.region} {params.annot_mp} {params.extra} 2> {output.logfile} |
            bcftools call -o {output.vcf} --multiallelic-caller --variants-only {params.ploidy} {params.annot_call} {params.groups} 2> {log}
        """

rule sort_genotypes:
    input:
        bcf = temp("call/{part}.vcf")
    output:
        bcf = temp("sort/{part}.bcf"),
        idx = temp("sort/{part}.bcf.csi")
    log:
        "logs/sort/{part}.sort.log"
    container:
        None
    shell:
        "bcftools sort --output {output.bcf} --write-index {input.bcf} 2> {log}"

rule concat_list:
    input:
        bcfs = collect("sort/{part}.bcf", part = intervals),
    output:
        "logs/bcf.files"
    run:
        with open(output[0], "w") as fout:
            for bcf in input.bcfs:
                _ = fout.write(f"{bcf}\n")  

rule concat_logs:
    input:
        collect("logs/mpileup/{part}.mpileup.log", part = intervals)
    output:
        "logs/mpileup.log"
    run:
        with open(output[0], "w") as fout:
            for file in input:
                interval = os.path.basename(file).replace(".mpileup.log", "")
                with open(file, "r") as fin:
                    for line in fin:
                        fout.write(f"{interval}\t{line}")

rule concat_variants:
    input:
        collect("sort/{part}.{ext}", part = intervals, ext = ["bcf", "bcf.csi"]),
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

rule indel_realign:
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
    container:
        None
    shell:
        "bcftools norm --threads {threads} -m -both -d both --write-index -Ob -o {output.bcf} -f {input.genome} {input.bcf} 2> {log}"

rule general_stats:
    input:
        genome  = workflow_geno,
        bcf     = "variants.{type}.bcf",
        idx     = "variants.{type}.bcf.csi"
    output:
        "reports/data/variants.{type}.stats"
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
        vcf = collect("variants.{file}.bcf", file = ["raw", "normalized"]),
        agg_log = "logs/mpileup.log",
        reports = collect("reports/variants.{file}.html", file = ["raw", "normalized"]) if not skip_reports else []
    params:
        ploidy = f"--ploidy {ploidy}",
        populations = f"--populations {groupings}" if groupings else "--populations -"
    run:
        summary = ["The harpy snp freebayes workflow ran using these parameters:"]
        summary.append(f"The provided reference genome: {bn}")
        summary.append(f"Genomic positions for which variants were called: {region_input}")
        mpileup = "The mpileup parameters:\n"
        mpileup += f"\tbcftools mpileup --fasta-ref REFERENCE --region REGION --bam-list BAMS --annotate AD --output-type b {mp_extra}"
        summary.append(mpileup)
        bcfcall = "The bcftools call parameters:\n"
        bcfcall += f"\tbcftools call --multiallelic-caller {params} --variants-only --output-type b |\n"
        bcfcall += "\tbcftools sort -"
        summary.append(bcfcall)
        merged = "The variants identified in the intervals were merged into the final variant file using:\n"
        merged += "\tbcftools concat -f bcf.files -a --remove-duplicates"
        summary.append(merged)
        normalize = "The variants were normalized using:\n"
        normalize += "\tbcftools norm -m -both -d both"
        summary.append(normalize)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['snakemake']['relative']}"
        summary.append(sm)
        with open("workflow/snp.mpileup.summary", "w") as f:
            f.write("\n\n".join(summary))