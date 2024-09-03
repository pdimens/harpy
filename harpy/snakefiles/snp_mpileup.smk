containerized: "docker://pdimens/harpy:latest"

import os
import sys
import multiprocessing
import logging as pylogging
from pathlib import Path

envdir      = os.getcwd() + "/.harpy_envs"
ploidy 		= config["ploidy"]
mp_extra 	= config.get("extra", "")
regiontype  = config["regiontype"]
windowsize  = config.get("windowsize", None)
outdir      = config["output_directory"]
skipreports = config["skip_reports"]
snakemake_log = config["snakemake_log"]
bamlist     = config["inputs"]["alignments"]
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
        while True:
            line = reg_in.readline()
            if not line:
                break
            cont,startpos,endpos = line.split()
            intervals.add(f"{cont}:{startpos}-{endpos}")
    regions = dict(zip(intervals, intervals))

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

onstart:
    extra_logfile_handler = pylogging.FileHandler(snakemake_log)
    logger.logger.addHandler(extra_logfile_handler)

def sam_index(infile):
    """Use Samtools to index an input file, adding .bai to the end of the name"""
    if not os.path.exists(f"{infile}.bai"):
        subprocess.run(f"samtools index {infile} {infile}.bai".split())

rule setup_genome:
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

rule faidx_genome:
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
        outdir + "/logs/samples.files"
    run:
        with open(output[0], "w") as fout:
            for bamfile in input.bam:
                _ = fout.write(bamfile + "\n")

rule mpileup:
    input:
        bamlist = outdir + "/logs/samples.files",
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
        outdir + "/variants.normalized.bcf"
    log:
        outdir + "/logs/variants.normalized.log"
    threads:
        workflow.cores
    container:
        None
    shell:
        "bcftools norm --threads {threads} -m -both -d both --write-index -Ob -o {output} -f {input.genome} {input.vcf} 2> {log}"

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
        outdir + "/reports/data/variants.{type}.stats"
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
        agg_log = outdir + "/logs/mpileup.log",
        reports = collect(outdir + "/reports/variants.{file}.html", file = ["raw", "normalized"]) if not skipreports else []
    params:
        ploidy = f"--ploidy {ploidy}",
        populations = f"--populations {groupings}" if groupings else "--populations -"
    run:
        import glob
        for logfile in glob.glob(f"{outdir}/logs/**/*", recursive = True):
            if os.path.isfile(logfile) and os.path.getsize(logfile) == 0:
                os.remove(logfile)
        for logfile in glob.glob(f"{outdir}/logs/**/*", recursive = True):
            if os.path.isdir(logfile) and not os.listdir(logfile):
                os.rmdir(logfile)
        with open(outdir + "/workflow/snp.mpileup.summary", "w") as f:
            _ = f.write("The harpy snp mpileup workflow ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            if windowsize:
                _ = f.write(f"Size of intervals to split genome for variant calling: {windowsize}\n")
            else:
                _ = f.write(f"Genomic positions for which variants were called: {regioninput}\n")
            _ = f.write("The mpileup parameters:\n")
            _ = f.write("    bcftools mpileup --fasta-ref GENOME --region REGION --bam-list BAMS --annotate AD --output-type b" + mp_extra + "\n")
            _ = f.write("The bcftools call parameters:\n")
            _ = f.write("    bcftools call --multiallelic-caller " + " ".join(params) + " --variants-only --output-type b | bcftools sort -\n")
            _ = f.write("The variants identified in the intervals were merged into the final variant file using:\n")
            _ = f.write("    bcftools concat -f bcf.files -a --remove-duplicates\n")
            _ = f.write("The variants were normalized using:\n")
            _ = f.write("    bcftools norm -m -both -d both\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")