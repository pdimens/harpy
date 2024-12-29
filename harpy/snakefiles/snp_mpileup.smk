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
mp_extra 	= config.get("extra", "")
regiontype  = config["region_type"]
windowsize  = config.get("windowsize", None)
skip_reports = config["reports"]["skip"]
bamlist     = config["inputs"]["alignments"]
bamdict     = dict(zip(bamlist, bamlist))
genomefile 	= config["inputs"]["genome"]
bn          = os.path.basename(genomefile)
genome_zip  = True if bn.lower().endswith(".gz") else False
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
        f"Genome/{bn}"
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
        f"Genome/{bn}"
    output: 
        fai = f"Genome/{bn}.fai",
        gzi = f"Genome/{bn}.gzi" if genome_zip else []
    log:
        f"Genome/{bn}.faidx.log"
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
        outdir + "/workflow/sample.groups"
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
        temp(outdir + "/logs/samples.files")
    run:
        with open(output[0], "w") as fout:
            for bamfile in input.bam:
                _ = fout.write(bamfile + "\n")

rule mpileup:
    input:
        bamlist = outdir + "/logs/samples.files",
        bam = bamlist,
        bai = collect("{bam}.bai", bam = bamlist),
        genome  = f"Genome/{bn}",
        genome_fai = f"Genome/{bn}.fai"
    output: 
        bcf = pipe(outdir + "/mpileup/{part}.mp.bcf"),
        logfile = temp(outdir + "/logs/{part}.mpileup.log")
    params:
        region = lambda wc: "-r " + regions[wc.part],
        extra = mp_extra
    container:
        None
    shell:
        "bcftools mpileup --fasta-ref {input.genome} --bam-list {input.bamlist} --annotate AD --output-type b {params} > {output.bcf} 2> {output.logfile}"

rule call_genotypes:
    input:
        groupfile = outdir + "/workflow/sample.groups" if groupings else [],
        bcf = outdir + "/mpileup/{part}.mp.bcf"
    output:
        bcf = temp(outdir + "/call/{part}.bcf"),
        idx = temp(outdir + "/call/{part}.bcf.csi")
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
        bcfs = collect(outdir + "/call/{part}.bcf", part = intervals),
    output:
        outdir + "/logs/bcf.files"
    run:
        with open(output[0], "w") as fout:
            for bcf in input.bcfs:
                _ = fout.write(f"{bcf}\n")  

rule concat_logs:
    input:
        collect(outdir + "/logs/{part}.mpileup.log", part = intervals)
    output:
        outdir + "/logs/mpileup.log"
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
        vcfs     = collect(outdir + "/call/{part}.{ext}", part = intervals, ext = ["bcf", "bcf.csi"]),
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
        bcf     = outdir + "/variants.{type}.bcf",
        idx     = outdir + "/variants.{type}.bcf.csi"
    output:
        outdir + "/reports/data/variants.{type}.stats"
    container:
        None
    shell:
        """
        bcftools stats -s "-" --fasta-ref {input.genome} {input.bcf} > {output} 2> /dev/null
        """

rule variant_report:
    input: 
        data = outdir + "/reports/data/variants.{type}.stats",
        qmd = f"{outdir}/workflow/report/bcftools_stats.qmd"
    output:
        report = outdir + "/reports/variants.{type}.html",
        qmd = temp(outdir + "/reports/variants.{type}.qmd")
    params:
        lambda wc: "variants." + wc.get("type")
    log:
        outdir + "/logs/variants.{type}.report.log"
    conda:
        f"{envdir}/r.yaml"
    shell:
        """
        cp {input.qmd} {output.qmd}
        INPATH=$(realpath {input.data})
        quarto render {output.qmd} -P infile:$INPATH inname:{params} 2> {log}
        """

rule workflow_summary:
    default_target: True
    input:
        vcf = collect(outdir + "/variants.{file}.bcf", file = ["raw"]),
        agg_log = outdir + "/logs/mpileup.log",
        reports = collect(outdir + "/reports/variants.{file}.html", file = ["raw", "normalized"]) if not skip_reports else []
    params:
        ploidy = f"--ploidy {ploidy}",
        populations = f"--populations {groupings}" if groupings else "--populations -"
    run:
        summary = ["The harpy snp freebayes workflow ran using these parameters:"]
        summary.append(f"The provided genome: {bn}")
        if windowsize:
            summary.append(f"Size of intervals to split genome for variant calling: {windowsize}")
        else:
            summary.append(f"Genomic positions for which variants were called: {regioninput}")
        mpileup = "The mpileup parameters:\n"
        mpileup += f"\tbcftools mpileup --fasta-ref GENOME --region REGION --bam-list BAMS --annotate AD --output-type b {mp_extra}"
        summary.append(mpileup)
        bcfcall = "The bcftools call parameters:\n"
        bcfcall += "\tbcftools call --multiallelic-caller {params} --variants-only --output-type b |\n"
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
        with open(outdir + "/workflow/snp.mpileup.summary", "w") as f:
            f.write("\n\n".join(summary))