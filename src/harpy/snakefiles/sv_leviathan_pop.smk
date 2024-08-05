containerized: "docker://pdimens/harpy:latest"

import os
import re
import sys
import logging as pylogging
from datetime import datetime
from rich import print as rprint
from rich.panel import Panel

envdir      = os.getcwd() + "/.harpy_envs"
genomefile 	= config["inputs"]["genome"]
bamlist     = config["inputs"]["alignments"]
groupfile 	= config["inputs"]["groupings"]
extra 		= config.get("extra", "") 
min_sv      = config["min_sv"]
min_bc      = config["min_barcodes"]
iterations  = config["iterations"]
outdir      = config["output_directory"]
skipreports = config["skip_reports"]
bn 			= os.path.basename(genomefile)
if bn.lower().endswith(".gz"):
    bn = bn[:-3]

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+",
    population = "[a-zA-Z0-9._-]+"
    
onstart:
    os.makedirs(f"{outdir}/logs/snakemake", exist_ok = True)
    dt_string = datetime.now().strftime("%d_%m_%Y-%H_%M_%S")
    extra_logfile_handler = pylogging.FileHandler(f"{outdir}/logs/snakemake/{dt_string}.snakelog")
    logger.logger.addHandler(extra_logfile_handler)

onerror:
    print("")
    rprint(
        Panel(
            f"The workflow has terminated due to an error. See the log file for more details:\n[bold]{outdir}/logs/snakemake/{dt_string}.snakelog[/bold]",
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

# create dictionary of population => filenames
## this makes it easier to set the snakemake rules/wildcards
def pop_manifest(groupingfile, filelist):
    d = {}
    with open(groupingfile) as f:
        for line in f:
            samp, pop = line.rstrip().split()
            if samp.lstrip().startswith("#"):
                continue
            r = re.compile(fr".*/({samp.lstrip()})\.(bam|sam)$", flags = re.IGNORECASE)
            sampl = list(filter(r.match, filelist))[0]
            if pop not in d.keys():
                d[pop] = [sampl]
            else:
                d[pop].append(sampl)
    return d

popdict = pop_manifest(groupfile, bamlist)
populations = popdict.keys()

rule copy_groupings:
    input:
        groupfile
    output:
        outdir + "/workflow/sample.groups"
    message:
        "Logging {input}"
    run:
        with open(input[0], "r") as infile, open(output[0], "w") as outfile:
            _ = [outfile.write(i) for i in infile.readlines() if not i.lstrip().startswith("#")]

rule merge_list:
    input:
        outdir + "/workflow/sample.groups"
    output:
        outdir + "/workflow/merge_samples/{population}.list"
    message:
        "Creating population file list: {wildcards.population}"
    run:
        with open(output[0], "w") as fout:
            for bamfile in popdict[wildcards.population]:
                _ = fout.write(bamfile + "\n")

rule merge_populations:
    input: 
        bamlist  = outdir + "/workflow/merge_samples/{population}.list",
        bamfiles = lambda wc: collect("{sample}", sample = popdict[wc.population]) 
    output:
        temp(outdir + "/workflow/input/{population}.unsort.bam")
    threads:
        1
    container:
        None
    message:
        "Merging alignments: {wildcards.population}"
    shell:
        "concatenate_bam.py -o {output} -b {input.bamlist}"

rule sort_merged:
    input:
        outdir + "/workflow/input/{population}.unsort.bam"
    output:
        bam = temp(outdir + "/workflow/input/{population}.bam"),
        bai = temp(outdir + "/workflow/input/{population}.bam.bai")
    log:
        outdir + "/logs/{population}.sort.log"
    resources:
        mem_mb = 2000
    threads:
        10
    container:
        None
    message:
        "Sorting alignments: {wildcards.population}"
    shell:
        "samtools sort -@ {threads} -O bam -l 0 -m {resources.mem_mb}M --write-index -o {output.bam}##idx##{output.bai} {input} 2> {log}"

rule index_barcode:
    input: 
        bam = outdir + "/workflow/input/{population}.bam",
        bai = outdir + "/workflow/input/{population}.bam.bai"
    output:
        temp(outdir + "/lrez_index/{population}.bci")
    benchmark:
        ".Benchmark/leviathan-pop/{population}.lrez"
    threads:
        max(10, workflow.cores)
    conda:
        f"{envdir}/sv.yaml"
    message:
        "Indexing barcodes: Population {wildcards.population}"
    shell:
        "LRez index bam -p -b {input.bam} -o {output} --threads {threads}"

rule genome_setup:
    input:
        genomefile
    output: 
        f"Genome/{bn}"
    container:
        None
    message: 
        "Creating {output}"
    shell: 
        "seqtk seq {input} > {output}"

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
    message:
        "Indexing {input}"
    conda:
        f"{envdir}/align.yaml"
    shell: 
        "bwa index {input} 2> {log}"

rule call_sv:
    input:
        bam    = outdir + "/workflow/input/{population}.bam",
        bai    = outdir + "/workflow/input/{population}.bam.bai",
        bc_idx = outdir + "/lrez_index/{population}.bci",
        genome = f"Genome/{bn}",
        genidx = multiext(f"Genome/{bn}", ".fai", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output:
        temp(outdir + "/vcf/{population}.vcf")
    log:  
        runlog     = outdir + "/logs/{population}.leviathan.log",
        candidates = outdir + "/logs/{population}.candidates"
    params:
        min_sv = f"-v {min_sv}",
        min_bc = f"-c {min_bc}",
        iters  = f"-B {iterations}",
        extra = extra
    threads:
        workflow.cores
    conda:
        f"{envdir}/sv.yaml"
    message:
        "Calling variants: {wildcards.population}"
    benchmark:
        ".Benchmark/leviathan-pop/{population}.variantcall"
    shell:
        "LEVIATHAN -b {input.bam} -i {input.bc_idx} {params} -g {input.genome} -o {output} -t {threads} --candidates {log.candidates} 2> {log.runlog}"

rule sort_bcf:
    input:
        outdir + "/vcf/{population}.vcf"
    output:
        outdir + "/vcf/{population}.bcf"
    params:
        "{wildcards.population}"
    container:
        None
    message:
        "Sorting and converting to BCF: Population {wildcards.population}"
    shell:        
        "bcftools sort -Ob --output {output} {input} 2> /dev/null"

rule sv_stats:
    input: 
        outdir + "/vcf/{population}.bcf"
    output:
        temp(outdir + "/reports/data/{population}.sv.stats")
    container:
        None
    message:
        "Getting stats: Population {input}"
    shell:
        """
        echo -e "population\\tcontig\\tposition_start\\tposition_end\\tlength\\ttype\\tn_barcodes\\tn_pairs" > {output}
        bcftools query -f '{wildcards.population}\\t%CHROM\\t%POS\\t%END\\t%SVLEN\\t%SVTYPE\\t%BARCODES\\t%PAIRS\\n' {input} >> {output}
        """

rule merge_variants:
    input:
        collect(outdir + "/reports/data/{population}.sv.stats", population = populations)
    output:
        outdir + "/inversions.bedpe",
        outdir + "/deletions.bedpe",
        outdir + "/duplications.bedpe",
        outdir + "/breakends.bedpe"
    message:
        "Aggregating the detected variants"
    run:
        with open(output[0], "w") as inversions, open(output[1], "w") as deletions, open(output[2], "w") as duplications, open(output[3], "w") as breakends:
            header = ["population","contig","position_start","position_end","length","type","n_barcodes","n_pairs"]
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

rule sv_report:
    input:	
        statsfile = outdir + "/reports/data/{population}.sv.stats",
        bcf       = outdir + "/vcf/{population}.bcf",
        faidx     = f"Genome/{bn}.fai"
    output:
        outdir + "/reports/{population}.sv.html"
    message:
        "Generating SV report: population {wildcards.population}"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/leviathan.Rmd"

rule sv_report_aggregate:
    input:	
        faidx      = f"Genome/{bn}.fai",
        statsfiles = collect(outdir + "/reports/data/{pop}.sv.stats", pop = populations)
    output:
        outdir + "/reports/leviathan.summary.html"
    message:
        "Generating SV report for all populations"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/leviathan_pop.Rmd"

rule workflow_summary:
    default_target: True
    input:
        vcf = collect(outdir + "/vcf/{pop}.bcf", pop = populations),
        bedpe_agg = collect(outdir + "/{sv}.bedpe", sv = ["inversions", "deletions","duplications", "breakends"]),
        reports = collect(outdir + "/reports/{pop}.sv.html", pop = populations) if not skipreports else [],
        agg_report = outdir + "/reports/leviathan.summary.html" if not skipreports else []
    message:
        "Summarizing the workflow: {output}"
    params:
        min_sv = f"-v {min_sv}",
        min_bc = f"-c {min_bc}",
        iters  = f"-B {iterations}",
        extra = extra
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