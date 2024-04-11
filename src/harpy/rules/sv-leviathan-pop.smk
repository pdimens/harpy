from rich import print as rprint
from rich.panel import Panel
import sys
import os

bam_dir 	= config["seq_directory"]
genomefile 	= config["genomefile"]
samplenames = config["samplenames"]
extra 		= config.get("extra", "") 
groupfile 	= config["groupings"]
min_sv      = config["min_sv"]
min_bc      = config["min_barcodes"]
outdir      = config["output_directory"]
skipreports = config["skipreports"]
bn 			= os.path.basename(genomefile)
genome_zip  = True if bn.lower().endswith(".gz") else False
if genome_zip:
    bn = bn[:-3]

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"
    
# create dictionary of population => filenames
## this makes it easier to set the snakemake rules/wildcards
## exits with an error if the groupfile has samples not in the bam folder
def pop_manifest(infile, dirn, sampnames):
    d = dict()
    with open(infile) as f:
        for line in f:
            samp, pop = line.rstrip().split()
            if samp.lstrip().startswith("#"):
                continue
            samp = f"{dirn}/{samp}.bam"
            if pop not in d.keys():
                d[pop] = [samp]
            else:
                d[pop].append(samp)
    return d

popdict = pop_manifest(groupfile, bam_dir, samplenames)
populations = popdict.keys()

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

rule copy_groupings:
    input:
        groupfile
    output:
        outdir + "/logs/sample.groups"
    message:
        "Logging {input}"
    run:
        with open(input[0], "r") as infile, open(output[0], "w") as outfile:
            _ = [outfile.write(i) for i in infile.readlines() if not i.lstrip().startswith("#")]

rule bam_list:
    input:
        outdir + "/logs/sample.groups"
    output:
        expand(outdir + "/input/{pop}.list", pop = populations)
    message:
        "Creating population file lists."
    run:
        for p in populations:
            bamlist = popdict[p]
            with open(f"{outdir}/input/{p}.list", "w") as fout:
                for bamfile in bamlist:
                    _ = fout.write(bamfile + "\n")

rule merge_populations:
    input: 
        bamlist  = outdir + "/input/{population}.list",
        bamfiles = lambda wc: expand("{sample}", sample = popdict[wc.population]) 
    output:
        bam = temp(outdir + "/input/{population}.bam"),
        bai = temp(outdir + "/input/{population}.bam.bai")
    threads:
        2
    message:
        "Merging alignments: Population {wildcards.population}"
    shell:
        "samtools merge -o {output.bam}##idx##{output.bai} --threads {threads} --write-index -b {input.bamlist}"

rule index_barcode:
    input: 
        bam = outdir + "/input/{population}.bam",
        bai = outdir + "/input/{population}.bam.bai"
    output:
        temp(outdir + "/lrezIndexed/{population}.bci")
    message:
        "Indexing barcodes: Population {wildcards.population}"
    benchmark:
        ".Benchmark/Variants/leviathan-pop/indexbc.{population}.txt"
    threads:
        4
    conda:
        os.getcwd() + "/.harpy_envs/variants.sv.yaml"
    shell:
        "LRez index bam -p -b {input.bam} -o {output} --threads {threads}"

rule genome_link:
    input:
        genomefile
    output: 
        f"Genome/{bn}"
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
            ln -sr {input} {output}
        fi
        """

rule genome_faidx:
    input: 
        f"Genome/{bn}"
    output: 
        f"Genome/{bn}.fai"
    log:
        f"Genome/{bn}.faidx.log"
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
    message:
        "Indexing {input}"
    conda:
        os.getcwd() + "/.harpy_envs/align.yaml"
    shell: 
        "bwa index {input} 2> {log}"

rule leviathan_variantcall:
    input:
        bam    = outdir + "/input/{population}.bam",
        bai    = outdir + "/input/{population}.bam.bai",
        bc_idx = outdir + "/lrezIndexed/{population}.bci",
        genome = f"Genome/{bn}",
        genidx = multiext(f"Genome/{bn}", ".fai", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output:
        pipe(outdir + "/{population}.vcf")
    log:  
        runlog     = outdir + "/logs/{population}.leviathan.log",
        candidates = outdir + "/logs/{population}.candidates"
    params:
        min_sv = f"-v {min_sv}",
        min_bc = f"-c {min_bc}",
        extra = extra
    threads:
        3
    conda:
        os.getcwd() + "/.harpy_envs/variants.sv.yaml"
    message:
        "Calling variants: Population {wildcards.population}"
    benchmark:
        ".Benchmark/Variants/leviathan-pop/variantcall.{population}.txt"
    shell:
        "LEVIATHAN -b {input.bam} -i {input.bc_idx} {params} -g {input.genome} -o {output} -t {threads} --candidates {log.candidates} 2> {log.runlog}"

rule sort_bcf:
    input:
        outdir + "/{population}.vcf"
    output:
        outdir + "/{population}.bcf"
    message:
        "Sorting and converting to BCF: Population {wildcards.population}"
    params:
        "{wildcards.population}"
    benchmark:
        ".Benchmark/Variants/leviathan-pop/sortbcf.{population}.txt"
    shell:        
        "bcftools sort -Ob --output {output} {input} 2> /dev/null"

rule sv_stats:
    input: 
        outdir + "/{population}.bcf"
    output:
        outdir + "/reports/reports/{population}.sv.stats"
    message:
        "Getting stats: Population {input}"
    benchmark:
        ".Benchmark/Variants/leviathan-pop/stats.{population}.txt"
    shell:
        """
        echo -e "population\\tcontig\\tposition_start\\tposition_end\\tlength\\ttype\\tn_barcodes\\tn_pairs" > {output}
        bcftools query -f '{wildcards.population}\\t%CHROM\\t%POS\\t%END\\t%SVLEN\\t%SVTYPE\\t%BARCODES\\t%PAIRS\\n' {input} >> {output}
        """

rule sv_report_bypop:
    input:	
        statsfile = outdir + "/reports/reports/{population}.sv.stats",
        bcf       = outdir + "/{population}.bcf",
        faidx     = f"Genome/{bn}.fai"
    output:
        outdir + "/reports/{population}.sv.html"
    message:
        "Generating SV report: population {wildcards.population}"
    conda:
        os.getcwd() + "/.harpy_envs/r-env.yaml"
    script:
        "report/Leviathan.Rmd"


rule sv_report:
    input:	
        faidx      = f"Genome/{bn}.fai",
        statsfiles = expand(outdir + "/reports/reports/{pop}.sv.stats", pop = populations)
    output:
        outdir + "/reports/leviathan.pop.summary.html"
    message:
        "Generating SV report for all populations"
    conda:
        os.getcwd() + "/.harpy_envs/r-env.yaml"
    script:
        "report/LeviathanPop.Rmd"

rule log_workflow:
    default_target: True
    input:
        vcf = expand(outdir + "/{pop}.bcf", pop = populations),
        reports = expand(outdir + "/reports/{pop}.sv.html", pop = populations),
        agg_report = outdir + "/reports/leviathan.pop.summary.html" if not skipreports else []
    output:
        outdir + "/workflow/sv.leviathan.summary"
    message:
        "Summarizing the workflow: {output}"
    params:
        min_sv = f"-v {min_sv}",
        min_bc = f"-c {min_bc}",
        extra = extra
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy variants sv module ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write(f"The directory with alignments: {bam_dir}\n")
            _ = f.write("The barcodes were indexed using:\n")
            _ = f.write("    LRez index bam -p -b INPUT\n")
            _ = f.write("Leviathan was called using:\n")
            _ = f.write(f"    LEVIATHAN -b INPUT -i INPUT.BCI -g GENOME {params}\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")