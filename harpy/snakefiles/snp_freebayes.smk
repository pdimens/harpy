containerized: "docker://pdimens/harpy:latest"

import os
import sys
import gzip
import multiprocessing
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

envdir      = os.getcwd() + "/.harpy_envs"
ploidy 		= config["ploidy"]
extra 	    = config.get("extra", "") 
regiontype  = config["regiontype"]
windowsize  = config.get("windowsize", None)
outdir      = config["output_directory"]
skipreports = config["skip_reports"]
bamlist     = config["inputs"]["alignments"]
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
if regiontype == "region":
    intervals = [regioninput]
    regions = {f"{regioninput}" : f"{regioninput}"}
else:
    with open(regioninput, "r") as reg_in:
        intervals = set()
        while True:
            line = reg_in.readline()
            if not line:
                break
            cont,startpos,endpos = line.split()
            intervals.add(f"{cont}:{startpos}-{endpos}")
    regions = dict(zip(intervals, intervals))

def sam_index(infile):
    """Use Samtools to index an input file, adding .bai to the end of the name"""
    if not os.path.exists(f"{infile}.bai"):
        subprocess.run(f"samtools index {infile} {infile}.bai".split())

rule preproc_groups:
    input:
        groupings
    output:
        outdir + "/workflow/sample.groups"
    run:
        with open(input[0], "r") as infile, open(output[0], "w") as outfile:
            _ = [outfile.write(i) for i in infile.readlines() if not i.lstrip().startswith("#")]

rule setup_genome:
    input:
        genomefile
    output: 
        f"Genome/{bn}"
    container:
        None
    shell: 
        "seqtk seq {input} > {output}"

rule faidx_genome:
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
        bamlist
    output:
        [f"{i}.bai" for i in bamlist]
    threads:
        workflow.cores
    run:
        with multiprocessing.Pool(processes=threads) as pool:
            pool.map(sam_index, input)

rule bam_list:
    input: 
        bam = bamlist,
        bai = [f"{i}.bai" for i in bamlist]
    output:
        outdir + "/workflow/samples.files"
    run:
        with open(output[0], "w") as fout:
            for bamfile in input.bam:
                _ = fout.write(bamfile + "\n")

rule call_variants:
    input:
        bam = bamlist,
        bai = [f"{i}.bai" for i in bamlist],
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
            bcftools sort - --output {out.bcf} --write-index 2> /dev/null
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
        outdir + "/reports/variants.{type}.stats",
    container:
        None
    shell:
        """
        bcftools stats -s "-" --fasta-ref {input.genome} {input.bcf} > {output} 2> /dev/null
        """

rule variant_report:
    input:
        outdir + "/reports/variants.{type}.stats"
    output:
        outdir + "/reports/variants.{type}.html"
    log:
        logfile = outdir + "/logs/variants.{type}.report.log"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/bcftools_stats.Rmd"

rule workflow_summary:
    default_target: True
    input:
        vcf = collect(outdir + "/variants.{file}.bcf", file = ["raw"]),
        reports = collect(outdir + "/reports/variants.{file}.html", file = ["raw"]) if not skipreports else []
    params:
        ploidy = f"-p {ploidy}",
        populations = f"--populations {groupings}" if groupings else '',
        extra = extra
    run:
        import glob
        for logfile in glob.glob(f"{outdir}/logs/**/*", recursive = True):
            if os.path.isfile(logfile) and os.path.getsize(logfile) == 0:
                os.remove(logfile)
        for logfile in glob.glob(f"{outdir}/logs/**/*", recursive = True):
            if os.path.isdir(logfile) and not os.listdir(logfile):
                os.rmdir(logfile)
        with open(outdir + "/workflow/snp.freebayes.summary", "w") as f:
            _ = f.write("The harpy snp freebayes workflow ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write(f"Size of intervals to split genome for variant calling: {windowsize}\n")
            _ = f.write("The freebayes parameters:\n")
            _ = f.write("    freebayes -f GENOME -L samples.list -r REGION " + " ".join(params) + " | bcftools sort -\n")
            _ = f.write("The variants identified in the intervals were merged into the final variant file using:\n")
            _ = f.write("    bcftools concat -f bcf.files -a --remove-duplicates\n")
            _ = f.write("The variants were normalized using:\n")
            _ = f.write("    bcftools norm -m -both -d both\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")