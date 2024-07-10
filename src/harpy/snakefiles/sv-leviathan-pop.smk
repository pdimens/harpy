containerized: "docker://pdimens/harpy:latest"

import os
import re
import sys
from rich import print as rprint
from rich.panel import Panel

envdir      = os.getcwd() + "/.harpy_envs"
genomefile 	= config["inputs"]["genome"]
bamlist     = config["inputs"]["alignments"]
groupfile 	= config["inputs"]["groupings"]
extra 		= config.get("extra", "") 
min_sv      = config["min_sv"]
min_bc      = config["min_barcodes"]
outdir      = config["output_directory"]
skipreports = config["skip_reports"]
bn 			= os.path.basename(genomefile)
if bn.lower().endswith(".gz"):
    bn = bn[:-3]

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+",
    population = "[a-zA-Z0-9._-]+"
    
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

rule bam_list:
    input:
        outdir + "/workflow/sample.groups"
    output:
        collect(outdir + "/workflow/{pop}.list", pop = populations)
    message:
        "Creating population file lists"
    run:
        for p in populations:
            alnlist = popdict[p]
            with open(f"{outdir}/workflow/{p}.list", "w") as fout:
                for bamfile in alnlist:
                    _ = fout.write(bamfile + "\n")

rule merge_populations:
    input: 
        bamlist  = outdir + "/workflow/{population}.list",
        bamfiles = lambda wc: collect("{sample}", sample = popdict[wc.population]) 
    output:
        bam = temp(outdir + "/workflow/input/{population}.bam"),
        bai = temp(outdir + "/workflow/input/{population}.bam.bai")
    threads:
        4
    container:
        None
    message:
        "Merging alignments: Population {wildcards.population}"
    shell:
        "samtools merge -o {output.bam}##idx##{output.bai} --threads {threads} --write-index -b {input.bamlist}"

rule index_barcode:
    input: 
        bam = outdir + "/workflow/input/{population}.bam",
        bai = outdir + "/workflow/input/{population}.bam.bai"
    output:
        temp(outdir + "/lrez_index/{population}.bci")
    benchmark:
        ".Benchmark/leviathan-pop/{population}.lrez"
    threads:
        4
    conda:
        f"{envdir}/sv.yaml"
    message:
        "Indexing barcodes: Population {wildcards.population}"
    shell:
        "LRez index bam -p -b {input.bam} -o {output} --threads {threads}"

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
        f"{envdir}/sv.yaml"
    message:
        "Calling variants: Population {wildcards.population}"
    benchmark:
        ".Benchmark/leviathan-pop/{population}.variantcall"
    shell:
        "LEVIATHAN -b {input.bam} -i {input.bc_idx} {params} -g {input.genome} -o {output} -t {threads} --candidates {log.candidates} 2> {log.runlog}"

rule sort_bcf:
    input:
        outdir + "/{population}.vcf"
    output:
        outdir + "/{population}.bcf"
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
        outdir + "/{population}.bcf"
    output:
        outdir + "/reports/data/{population}.sv.stats"
    container:
        None
    message:
        "Getting stats: Population {input}"
    shell:
        """
        echo -e "population\\tcontig\\tposition_start\\tposition_end\\tlength\\ttype\\tn_barcodes\\tn_pairs" > {output}
        bcftools query -f '{wildcards.population}\\t%CHROM\\t%POS\\t%END\\t%SVLEN\\t%SVTYPE\\t%BARCODES\\t%PAIRS\\n' {input} >> {output}
        """

rule sv_report_bypop:
    input:	
        statsfile = outdir + "/reports/data/{population}.sv.stats",
        bcf       = outdir + "/{population}.bcf",
        faidx     = f"Genome/{bn}.fai"
    output:
        outdir + "/reports/{population}.sv.html"
    message:
        "Generating SV report: population {wildcards.population}"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/Leviathan.Rmd"

rule sv_report:
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
        "report/LeviathanPop.Rmd"

rule log_workflow:
    default_target: True
    input:
        vcf = collect(outdir + "/{pop}.bcf", pop = populations),
        reports = collect(outdir + "/reports/{pop}.sv.html", pop = populations) if not skipreports else [],
        agg_report = outdir + "/reports/leviathan.summary.html" if not skipreports else []
    message:
        "Summarizing the workflow: {output}"
    params:
        min_sv = f"-v {min_sv}",
        min_bc = f"-c {min_bc}",
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