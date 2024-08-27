containerized: "docker://pdimens/harpy:latest"

import os
import re
import sys
import multiprocessing
import logging as pylogging
from pathlib import Path

envdir      = os.getcwd() + "/.harpy_envs"
genomefile  = config["inputs"]["genome"]
bamlist     = config["inputs"]["alignments"]
samplenames = {Path(i).stem for i in bamlist}
min_sv      = config["min_sv"]
min_bc      = config["min_barcodes"]
iterations  = config["iterations"]
extra       = config.get("extra", "") 
outdir      = config["output_directory"]
skipreports = config["skip_reports"]
snakemake_log = config["snakemake_log"]
bn          = os.path.basename(genomefile)
genome_zip  = True if bn.lower().endswith(".gz") else False
if genome_zip:
    bn = bn[:-3]

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

onstart:
    extra_logfile_handler = pylogging.FileHandler(snakemake_log)
    logger.logger.addHandler(extra_logfile_handler)

def sam_index(infile):
    """Use Samtools to index an input file, adding .bai to the end of the name"""
    if not os.path.exists(f"{infile}.bai"):
        subprocess.run(f"samtools index {infile} {infile}.bai".split())

def get_alignments(wildcards):
    """returns a list with the bam file for the sample based on wildcards.sample"""
    r = re.compile(fr".*/({wildcards.sample})\.(bam|sam)$", flags = re.IGNORECASE)
    aln = list(filter(r.match, bamlist))
    return aln[0]

def get_align_index(wildcards):
    """returns a list with the bai index file for the sample based on wildcards.sample"""
    r = re.compile(fr"(.*/{wildcards.sample})\.(bam|sam)$", flags = re.IGNORECASE)
    aln = list(filter(r.match, bamlist))
    return aln[0] + ".bai"

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

rule index_barcode:
    input: 
        bam = get_alignments,
        bai = get_align_index
    output:
        temp(outdir + "/lrezIndexed/{sample}.bci")
    benchmark:
        ".Benchmark/leviathan/{sample}.lrez"
    threads:
        max(10, workflow.cores)
    conda:
        f"{envdir}/sv.yaml"
    shell:
        "LRez index bam --threads {threads} -p -b {input.bam} -o {output}"

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

rule bwa_index_genome:
    input: 
        f"Genome/{bn}"
    output: 
        multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    log:
        f"Genome/{bn}.bwa.idx.log"
    conda:
        f"{envdir}/align.yaml"
    shell: 
        "bwa index {input} 2> {log}"

rule call_variants:
    input:
        bam    = get_alignments,
        bai    = get_align_index,
        bc_idx = outdir + "/lrezIndexed/{sample}.bci",
        genome = f"Genome/{bn}",
        genidx = multiext(f"Genome/{bn}", ".fai", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output:
        vcf = temp(outdir + "/vcf/{sample}.vcf"),
        candidates = outdir + "/logs/leviathan/{sample}.candidates"
    log:  
        runlog = outdir + "/logs/leviathan/{sample}.leviathan.log"
    params:
        min_sv = f"-v {min_sv}",
        min_bc = f"-c {min_bc}",
        iters  = f"-B {iterations}",
        extra = extra
    threads:
        workflow.cores - 1
    conda:
        f"{envdir}/sv.yaml"
    benchmark:
        ".Benchmark/leviathan/{sample}.variantcall"
    shell:
        "LEVIATHAN -b {input.bam} -i {input.bc_idx} {params} -g {input.genome} -o {output.vcf} -t {threads} --candidates {output.candidates} 2> {log.runlog}"

rule sort_variants:
    input:
        outdir + "/vcf/{sample}.vcf"
    output:
        outdir + "/vcf/{sample}.bcf"
    params:
        "{wildcards.sample}"
    container:
        None
    shell:        
        "bcftools sort -Ob --output {output} {input} 2> /dev/null"

rule variant_stats:
    input: 
        outdir + "/vcf/{sample}.bcf"
    output: 
        temp(outdir + "/reports/data/{sample}.sv.stats")
    container:
        None
    shell:
        """
        echo -e "sample\\tcontig\\tposition_start\\tposition_end\\tlength\\ttype\\tn_barcodes\\tn_pairs" > {output}
        bcftools query -f '{wildcards.sample}\\t%CHROM\\t%POS\\t%END\\t%SVLEN\\t%SVTYPE\\t%BARCODES\\t%PAIRS\\n' {input} >> {output}
        """

rule aggregate_variants:
    input:
        collect(outdir + "/reports/data/{sample}.sv.stats", sample = samplenames)
    output:
        outdir + "/inversions.bedpe",
        outdir + "/deletions.bedpe",
        outdir + "/duplications.bedpe",
        outdir + "/breakends.bedpe"
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

rule sample_reports:
    input:	
        faidx     = f"Genome/{bn}.fai",
        statsfile = outdir + "/reports/data/{sample}.sv.stats"
    output:	
        outdir + "/reports/{sample}.SV.html"
    log:
        logfile = outdir + "/logs/reports/{sample}.report.log"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/leviathan.Rmd"

rule workflow_summary:
    default_target: True
    input: 
        vcf = collect(outdir + "/vcf/{sample}.bcf", sample = samplenames),
        bedpe_agg = collect(outdir + "/{sv}.bedpe", sv = ["inversions", "deletions","duplications", "breakends"]),
        reports = collect(outdir + "/reports/{sample}.SV.html", sample = samplenames) if not skipreports else []
    params:
        min_sv = f"-v {min_sv}",
        min_bc = f"-c {min_bc}",
        iters  = f"-B {iterations}",
        extra = extra
    run:
        import glob
        for logfile in glob.glob(f"{outdir}/logs/**/*", recursive = True):
            if os.path.isfile(logfile) and os.path.getsize(logfile) == 0:
                os.remove(logfile)
        for logfile in glob.glob(f"{outdir}/logs/**/*", recursive = True):
            if os.path.isdir(logfile) and not os.listdir(logfile):
                os.rmdir(logfile)
        with open(outdir + "/workflow/sv.leviathan.summary", "w") as f:
            _ = f.write("The harpy sv leviathan workflow ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write("The barcodes were indexed using:\n")
            _ = f.write("    LRez index bam -p -b INPUT\n")
            _ = f.write("Leviathan was called using:\n")
            _ = f.write(f"    LEVIATHAN -b INPUT -i INPUT.BCI -g GENOME {params}\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")