containerized: "docker://pdimens/harpy:latest"

import os
import logging
from pathlib import Path

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake_log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

envdir      = os.path.join(os.getcwd(), "workflow", "envs")
ploidy 		= config["ploidy"]
mp_extra 	= config.get("extra", "")
regiontype  = config["region_type"]
windowsize  = config.get("windowsize", None)
skip_reports = config["reports"]["skip"]
bamlist     = config["inputs"]["alignments"]
bamdict     = dict(zip(bamlist, bamlist))
genomefile 	= config["inputs"]["reference"]
bn          = os.path.basename(genomefile)
workflow_geno = f"workflow/reference/{bn}"
genome_zip  = True if bn.lower().endswith(".gz") else False
workflow_geno_idx = f"{workflow_geno}.gzi" if genome_zip else f"{workflow_geno}.fai"
groupings 	= config["inputs"].get("groupings", [])
regioninput = config["inputs"]["regions"]
samplenames = {Path(i).stem for i in bamlist}

if regiontype == "region":
    intervals = [regioninput]
    regions = {f"{regioninput}" : f"{regioninput}"}
else:
    with open(regioninput, "r") as reg_in:
        intervals = set()
        for line in reg_in:
            cont,startpos,endpos = line.split()
            intervals.add(f"{cont}:{startpos}-{endpos}")
    regions = dict(zip(intervals, intervals))

rule process_genome:
    input:
        genomefile
    output: 
        workflow_geno
    container:
        None
    shell: 
        """
        if (file {input} | grep -q compressed ) ;then
            # is regular gzipped, needs to be BGzipped
            seqtk seq {input} | bgzip -c > {output}
        else
            cp -f {input} {output}
        fi
        """

rule index_genome:
    input: 
        workflow_geno
    output: 
        fai = f"{workflow_geno}.fai",
        gzi = f"{workflow_geno}.gzi" if genome_zip else []
    log:
        f"{workflow_geno}.faidx.log"
    params:
        genome_zip
    container:
        None
    shell: 
        """
        if [ "{params}" = "True" ]; then
            samtools faidx --gzi-idx {output.gzi} --fai-idx {output.fai} {input} 2> {log}
        else
            samtools faidx --fai-idx {output.fai} {input} 2> {log}
        fi
        """

rule preproc_groups:
    input:
        groupings
    output:
        "workflow/sample.groups"
    run:
        with open(input[0], "r") as infile, open(output[0], "w") as outfile:
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
        bam = bamlist,
        bai = collect("{bam}.bai", bam = bamlist)
    output:
        temp("logs/samples.files")
    run:
        with open(output[0], "w") as fout:
            for bamfile in input.bam:
                _ = fout.write(bamfile + "\n")

rule mpileup:
    input:
        bamlist = "logs/samples.files",
        bam = bamlist,
        bai = collect("{bam}.bai", bam = bamlist),
        genome  = workflow_geno,
        genome_fai = f"{workflow_geno}.fai"
    output: 
        bcf = pipe("mpileup/{part}.mp.bcf"),
        logfile = temp("logs/{part}.mpileup.log")
    params:
        region = lambda wc: "-r " + regions[wc.part],
        extra = mp_extra
    container:
        None
    shell:
        "bcftools mpileup --fasta-ref {input.genome} --bam-list {input.bamlist} --annotate AD --output-type b {params} > {output.bcf} 2> {output.logfile}"

rule call_genotypes:
    input:
        groupfile = "workflow/sample.groups" if groupings else [],
        bcf = "mpileup/{part}.mp.bcf"
    output:
        bcf = temp("call/{part}.bcf"),
        idx = temp("call/{part}.bcf.csi")
    params: 
        f"--ploidy {ploidy}",
        "--group-samples" if groupings else "--group-samples -"
    threads:
        2
    container:
        None
    shell:
        """
        bcftools call --multiallelic-caller --variants-only --output-type b {params} {input} |
            bcftools sort - --output {output.bcf} --write-index 2> /dev/null
        """

rule concat_list:
    input:
        bcfs = collect("call/{part}.bcf", part = intervals),
    output:
        "logs/bcf.files"
    run:
        with open(output[0], "w") as fout:
            for bcf in input.bcfs:
                _ = fout.write(f"{bcf}\n")  

rule concat_logs:
    input:
        collect("logs/{part}.mpileup.log", part = intervals)
    output:
        "logs/mpileup.log"
    run:
        with open(output[0], "w") as fout:
            for file in input:
                fin = open(file, "r")
                interval = os.path.basename(file).replace(".mpileup.log", "")
                for line in fin.readlines():
                    fout.write(f"{interval}\t{line}")
                fin.close()

rule concat_variants:
    input:
        vcfs     = collect("call/{part}.{ext}", part = intervals, ext = ["bcf", "bcf.csi"]),
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

rule report_config:
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
        f"{envdir}/r.yaml"
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
        if windowsize:
            summary.append(f"Size of intervals to split genome for variant calling: {windowsize}")
        else:
            summary.append(f"Genomic positions for which variants were called: {regioninput}")
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
        sm += f"\t{config['workflow_call']}"
        summary.append(sm)
        with open("workflow/snp.mpileup.summary", "w") as f:
            f.write("\n\n".join(summary))