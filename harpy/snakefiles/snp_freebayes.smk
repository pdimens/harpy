containerized: "docker://pdimens/harpy:latest"

import os
import logging
from pathlib import Path

onstart:
    logger.logger.addHandler(logging.FileHandler(config["snakemake_log"]))
onsuccess:
    os.remove(logger.logfile)
onerror:
    os.remove(logger.logfile)
wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

outdir      = config["output_directory"]
envdir      = os.path.join(os.getcwd(), outdir, "workflow", "envs")
ploidy 		= config["ploidy"]
extra 	    = config.get("extra", "") 
regiontype  = config["region_type"]
windowsize  = config.get("windowsize", None)
skip_reports = config["reports"]["skip"]
bamlist     = config["inputs"]["alignments"]
bamdict     = dict(zip(bamlist, bamlist))
genomefile 	= config["inputs"]["genome"]
bn          = os.path.basename(genomefile)
if bn.lower().endswith(".gz"):
    genome_zip  = True
    bn = bn[:-3]
else:
    genome_zip  = False
groupings 	= config["inputs"].get("groupings", [])
regioninput = config["inputs"]["regions"]
samplenames = {Path(i).stem for i in bamlist}
sampldict = dict(zip(bamlist, samplenames))
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

rule preproc_groups:
    input:
        groupings
    output:
        outdir + "/workflow/sample.groups"
    run:
        with open(input[0], "r") as infile, open(output[0], "w") as outfile:
            _ = [outfile.write(i) for i in infile.readlines() if not i.lstrip().startswith("#")]

rule process_genome:
    input:
        genomefile
    output: 
        f"Genome/{bn}"
    container:
        None
    shell: 
        "seqtk seq {input} > {output}"

rule index_genome:
    input: 
        f"Genome/{bn}"
    output: 
        f"Genome/{bn}.fai"
    log:
        f"Genome/{bn}.faidx.log"
    container:
        None
    shell:
        "samtools faidx --fai-idx {output} {input} 2> {log}"

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
        outdir + "/workflow/samples.files"
    run:
        with open(output[0], "w") as fout:
            for bamfile in input.bam:
                _ = fout.write(bamfile + "\n")

rule call_variants:
    input:
        bam = bamlist,
        bai = collect("{bam}.bai", bam = bamlist),
        groupfile = outdir + "/workflow/sample.groups" if groupings else [],
        ref     = f"Genome/{bn}",
        ref_idx = f"Genome/{bn}.fai",
        samples = outdir + "/workflow/samples.files"
    output:
        bcf = temp(outdir + "/regions/{part}.bcf"),
        idx = temp(outdir + "/regions/{part}.bcf.csi")
    log:
        outdir + "/logs/{part}.freebayes.log"
    params:
        region = lambda wc: "-r " + regions[wc.part],
        ploidy = f"-p {ploidy}",
        populations = f"--populations {outdir}/workflow/sample.groups" if groupings else "",
        extra = extra
    threads:
        2
    conda:
        f"{envdir}/variants.yaml"
    shell:
        """
        freebayes -f {input.ref} -L {input.samples} {params} 2> {log} |
            bcftools sort - --output {output.bcf} --write-index 2> /dev/null
        """

rule concat_list:
    input:
        bcfs = collect(outdir + "/regions/{part}.bcf", part = intervals),
    output:
        outdir + "/logs/bcf.files"
    run:
        with open(output[0], "w") as fout:
            for bcf in input.bcfs:
                _ = fout.write(f"{bcf}\n")

rule concat_variants:
    input:
        bcfs = collect(outdir + "/regions/{part}.{ext}", part = intervals, ext = ["bcf", "bcf.csi"]),
        filelist = outdir + "/logs/bcf.files"
    output:
        temp(outdir + "/variants.raw.unsort.bcf")
    log:
        outdir + "/logs/concat.log"
    threads:
        workflow.cores
    container:
        None
    shell:  
        "bcftools concat -f {input.filelist} --threads {threads} --naive -Ob -o {output} 2> {log}"

rule sort_variants:
    input:
        outdir + "/variants.raw.unsort.bcf"
    output:
        bcf = outdir + "/variants.raw.bcf",
        csi = outdir + "/variants.raw.bcf.csi"
    container:
        None
    shell:
        "bcftools sort --write-index -Ob -o {output.bcf} {input} 2> /dev/null"

rule indel_realign:
    input:
        genome  = f"Genome/{bn}",
        bcf     = outdir + "/variants.raw.bcf",
        idx     = outdir + "/variants.raw.bcf.csi"
    output:
        bcf = outdir + "/variants.normalized.bcf",
        idx = outdir + "/variants.normalized.bcf.csi"
    log:
        outdir + "/logs/variants.normalized.log"
    threads:
        workflow.cores
    container:
        None
    shell:
        "bcftools norm --threads {threads} -m -both -d both --write-index -Ob -o {output.bcf} -f {input.genome} {input.bcf} 2> {log}"

rule general_stats:
    input:
        genome  = f"Genome/{bn}",
        ref_idx = f"Genome/{bn}.fai",
        bcf     = outdir + "/variants.{type}.bcf",
        idx     = outdir + "/variants.{type}.bcf.csi"
    output:
        outdir + "/reports/data/variants.{type}.stats",
    container:
        None
    shell:
        """
        bcftools stats -s "-" --fasta-ref {input.genome} {input.bcf} > {output} 2> /dev/null
        """

rule report_config:
    input:
        f"{outdir}/workflow/report/_quarto.yml"
    output:
        f"{outdir}/reports/_quarto.yml"
    shell:
        "cp {input} {output}"

rule variant_report:
    input: 
        data = outdir + "/reports/data/variants.{type}.stats",
        qmd  = f"{outdir}/workflow/report/bcftools_stats.qmd",
        yml  = f"{outdir}/reports/_quarto.yml"
    output:
        report = outdir + "/reports/variants.{type}.html",
        qmd = temp(outdir + "/reports/variants.{type}.qmd")
    params:
        lambda wc: "-P vcf:variants." + wc.get("type")
    log:
        outdir + "/logs/variants.{type}.report.log"
    conda:
        f"{envdir}/r.yaml"
    shell:
        """
        cp {input.qmd} {output.qmd}
        INPATH=$(realpath {input.data})
        quarto render {output.qmd} --log {log} --quiet -P infile:$INPATH -P vcf:{params}
        """

rule workflow_summary:
    default_target: True
    input:
        vcf = collect(outdir + "/variants.{file}.bcf", file = ["raw"]),
        reports = collect(outdir + "/reports/variants.{file}.html", file = ["raw"]) if not skip_reports else []
    params:
        ploidy = f"-p {ploidy}",
        populations = f"--populations {groupings}" if groupings else '',
        extra = extra
    run:
        summary = ["The harpy snp freebayes workflow ran using these parameters:"]
        summary.append(f"The provided genome: {bn}")
        if windowsize:
            summary.append(f"Size of intervals to split genome for variant calling: {windowsize}")
        else:
            summary.append(f"Genomic positions for which variants were called: {regioninput}")
        varcall = "The freebayes parameters:\n"
        varcall += f"\tfreebayes -f GENOME -L samples.list -r REGION {params} |\n"
        varcall += f"\tbcftools sort -"
        summary.append(varcall)
        merged = "The variants identified in the intervals were merged into the final variant file using:\n"
        merged += "\tbcftools concat -f bcf.files -a --remove-duplicates"
        summary.append(merged)
        normalize = "The variants were normalized using:\n"
        normalize += "\tbcftools norm -m -both -d both"
        summary.append(normalize)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['workflow_call']}"
        summary.append(sm)
        with open(outdir + "/workflow/snp.freebayes.summary", "w") as f:
            f.write("\n\n".join(summary))
