containerized: "docker://pdimens/harpy:latest"

import os
import re
import sys
import multiprocessing
from pathlib import Path
from rich import print as rprint
from rich.panel import Panel

envdir      = os.getcwd() + "/.harpy_envs"
genomefile  = config["inputs"]["genome"]
bamlist     = config["inputs"]["alignments"]
samplenames = {Path(i).stem for i in bamlist}
min_sv      = config["min_sv"]
min_bc      = config["min_barcodes"]
extra       = config.get("extra", "") 
outdir      = config["output_directory"]
skipreports = config["skipreports"]
bn          = os.path.basename(genomefile)
genome_zip  = True if bn.lower().endswith(".gz") else False
if genome_zip:
    bn = bn[:-3]

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

onerror:
    print("")
    rprint(
        Panel(
            f"The workflow has terminated due to an error. See the log file below for more details.",
            title = "[bold]harpy sv leviathan",
            title_align = "left",
            border_style = "red"
            ),
        file = sys.stderr
    )

onsuccess:
    print("")
    rprint(
        Panel(
            f"The workflow has finished successfully! Find the results in [bold]{outdir}/[/bold]",
            title = "[bold]harpy sv leviathan",
            title_align = "left",
            border_style = "green"
            ),
        file = sys.stderr
    )

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
    message:
        "Indexing alignment files"
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
        4
    conda:
        f"{envdir}/sv.yaml"
    message:
        "Indexing barcodes: {wildcards.sample}"
    shell:
        "LRez index bam --threads {threads} -p -b {input.bam} -o {output}"

rule genome_link:
    input:
        genomefile
    output: 
        f"Genome/{bn}"
    container:
        None
    message: 
        "Creating {output}"
    shell: 
        """
        if (file {input} | grep -q compressed ) ;then
            # is regular gzipped, needs to be decompressed
            gzip -dc {input} > {output}
        elif (file {input} | grep -q BGZF ); then
            # is bgzipped, decompress
            gzip -dc {input} > {output}
        else
            # isn't compressed, just linked
            cp -f {input} {output}
        fi
        """

rule genome_faidx:
    input: 
        f"Genome/{bn}"
    output: 
        f"Genome/{bn}.fai"
    log:
        f"Genome/{bn}.faidx.log"
    container:
        None
    message:
        "Indexing {input}"
    shell:
        "samtools faidx --fai-idx {output} {input} 2> {log}"

rule index_bwa_genome:
    input: 
        f"Genome/{bn}"
    output: 
        multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    log:
        f"Genome/{bn}.idx.log"
    conda:
        f"{envdir}/align.yaml"
    message:
        "Indexing {input}"
    shell: 
        "bwa index {input} 2> {log}"

rule call_sv:
    input:
        bam    = get_alignments,
        bai    = get_align_index,
        bc_idx = outdir + "/lrezIndexed/{sample}.bci",
        genome = f"Genome/{bn}",
        genidx = multiext(f"Genome/{bn}", ".fai", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output:
        pipe(outdir + "/{sample}.vcf")
    log:  
        runlog     = outdir + "/logs/{sample}.leviathan.log",
        candidates = outdir + "/logs/{sample}.candidates"
    params:
        min_sv = f"-v {min_sv}",
        min_bc = f"-c {min_bc}",
        extra = extra
    threads:
        3
    conda:
        f"{envdir}/sv.yaml"
    benchmark:
        ".Benchmark/leviathan/{sample}.variantcall"
    message:
        "Calling variants: {wildcards.sample}"
    shell:
        "LEVIATHAN -b {input.bam} -i {input.bc_idx} {params} -g {input.genome} -o {output} -t {threads} --candidates {log.candidates} 2> {log.runlog}"

rule sort_bcf:
    input:
        outdir + "/{sample}.vcf"
    output:
        outdir + "/{sample}.bcf"
    params:
        "{wildcards.sample}"
    container:
        None
    message:
        "Sorting and converting to BCF: {wildcards.sample}"
    shell:        
        "bcftools sort -Ob --output {output} {input} 2> /dev/null"

rule sv_stats:
    input: 
        outdir + "/{sample}.bcf"
    output: 
        outdir + "/reports/data/{sample}.sv.stats"
    container:
        None
    message:
        "Getting SV stats for {wildcards.sample}"
    shell:
        """
        echo -e "sample\\tcontig\\tposition_start\\tposition_end\\tlength\\ttype\\tn_barcodes\\tn_pairs" > {output}
        bcftools query -f '{wildcards.sample}\\t%CHROM\\t%POS\\t%END\\t%SVLEN\\t%SVTYPE\\t%BARCODES\\t%PAIRS\\n' {input} >> {output}
        """

rule sv_report:
    input:	
        faidx     = f"Genome/{bn}.fai",
        statsfile = outdir + "/reports/data/{sample}.sv.stats"
    output:	
        outdir + "/reports/{sample}.SV.html"
    conda:
        f"{envdir}/r.yaml"
    message:
        "Generating SV report: {wildcards.sample}"
    script:
        "report/Leviathan.Rmd"


rule log_workflow:
    default_target: True
    input: 
        vcf = collect(outdir + "/{sample}.bcf", sample = samplenames),
        reports = collect(outdir + "/reports/{sample}.SV.html", sample = samplenames) if not skipreports else []
    params:
        min_sv = f"-v {min_sv}",
        min_bc = f"-c {min_bc}",
        extra = extra
    message:
        "Summarizing the workflow: {output}"
    run:
        with open(outdir + "/workflow/sv.leviathan.summary", "w") as f:
            _ = f.write("The harpy sv leviathan workflow ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write("The barcodes were indexed using:\n")
            _ = f.write("    LRez index bam -p -b INPUT\n")
            _ = f.write("Leviathan was called using:\n")
            _ = f.write(f"    LEVIATHAN -b INPUT -i INPUT.BCI -g GENOME {params}\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")