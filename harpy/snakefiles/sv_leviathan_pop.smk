containerized: "docker://pdimens/harpy:latest"

import os
import re
import logging

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake"]["log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+",
    population = r"[a-zA-Z0-9._-]+"

genomefile 	= config["inputs"]["reference"]
bamlist     = config["inputs"]["alignments"]
groupfile 	= config["inputs"]["groupings"]
extra 		= config.get("extra", "") 
min_size      = config["min_size"]
min_bc      = config["min_barcodes"]
iterations  = config["iterations"]
small_thresh = config["variant_thresholds"]["small"]
medium_thresh = config["variant_thresholds"]["medium"]
large_thresh = config["variant_thresholds"]["large"]
duplcates_thresh = config["variant_thresholds"]["duplicates"]
skip_reports = config["reports"]["skip"]
plot_contigs = config["reports"]["plot_contigs"]
bn 			= os.path.basename(genomefile)
if bn.lower().endswith(".gz"):
    workflow_geno = f"workflow/reference/{bn[:-3]}"
else:
    workflow_geno = f"workflow/reference/{bn}"
    
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

rule preprocess_groups:
    input:
        grp = groupfile
    output:
        grp = "workflow/sample.groups"
    run:
        with open(input.grp, "r") as infile, open(output.grp, "w") as outfile:
            _ = [outfile.write(i) for i in infile.readlines() if not i.lstrip().startswith("#")]

rule concat_list:
    input:
        "workflow/sample.groups"
    output:
        "workflow/merge_samples/{population}.list"
    run:
        with open(output[0], "w") as fout:
            for bamfile in popdict[wildcards.population]:
                _ = fout.write(bamfile + "\n")

rule concat_groups:
    input: 
        bamlist  = "workflow/merge_samples/{population}.list",
        bamfiles = lambda wc: collect("{sample}", sample = popdict[wc.population]) 
    output:
        temp("workflow/input/{population}.unsort.bam")
    log:
        "logs/concat_groups/{population}.concat.log"
    threads:
        1
    container:
        None
    shell:
        "concatenate_bam.py --bx -b {input.bamlist} > {output} 2> {log}"

rule sort_groups:
    input:
        "workflow/input/{population}.unsort.bam"
    output:
        bam = temp("workflow/input/{population}.bam"),
        bai = temp("workflow/input/{population}.bam.bai")
    log:
        "logs/samtools_sort/{population}.sort.log"
    resources:
        mem_mb = 2000
    threads:
        10
    container:
        None
    shell:
        "samtools sort -@ {threads} -O bam -l 0 -m {resources.mem_mb}M --write-index -o {output.bam}##idx##{output.bai} {input} 2> {log}"

rule index_barcode:
    input: 
        bam = "workflow/input/{population}.bam",
        bai = "workflow/input/{population}.bam.bai"
    output:
        temp("lrez_index/{population}.bci")
    threads:
        min(5, workflow.cores)
    conda:
        "envs/variants.yaml"
    shell:
        "LRez index bam -p -b {input.bam} -o {output} --threads {threads}"

rule preprocess_reference:
    input:
        genomefile
    output: 
        multiext(workflow_geno, ".ann", ".bwt", ".pac", ".sa", ".amb"),
        geno = workflow_geno,
        fai  = f"{workflow_geno}.fai"
    log:
        f"{workflow_geno}.preprocess.log"
    conda:
        "envs/align.yaml"
    shell: 
        """
        seqtk seq {input} > {output.geno}
        samtools faidx --fai-idx {output.fai} {output.geno} 2> {log}
        bwa index {output.geno} 2>> {log}
        """

rule call_variants:
    input:
        bam    = "workflow/input/{population}.bam",
        bai    = "workflow/input/{population}.bam.bai",
        bc_idx = "lrez_index/{population}.bci",
        genome = workflow_geno,
        genidx = multiext(workflow_geno, ".fai", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output:
        vcf = temp("vcf/{population}.vcf"),
        candidates = "logs/leviathan/{population}.candidates"
    log:  
        runlog = "logs/leviathan/{population}.leviathan.log",
    params:
        min_size = f"-v {min_size}",
        min_bc = f"-c {min_bc}",
        iters  = f"-B {iterations}",
        small  = f"-s {small_thresh}",
        medium  = f"-m {medium_thresh}",
        large  = f"-l {large_thresh}",
        dupes  = f"-d {duplcates_thresh}",
        extra = extra
    threads:
        workflow.cores - 1
    conda:
        "envs/variants.yaml"
    shell:
        "LEVIATHAN -b {input.bam} -i {input.bc_idx} {params} -g {input.genome} -o {output.vcf} -t {threads} --candidates {output.candidates} 2> {log.runlog}"

rule sort_variants:
    priority: 100
    input:
        "vcf/{population}.vcf"
    output:
        "vcf/{population}.bcf"
    params:
        lambda wc: wc.population
    container:
        None
    shell:        
        "bcftools sort -Ob --output {output} {input} 2> /dev/null"

rule variant_stats:
    input: 
        "vcf/{population}.bcf"
    output:
        temp("reports/data/{population}.sv.stats")
    container:
        None
    shell:
        """
        echo -e "population\\tcontig\\tposition_start\\tposition_end\\tlength\\ttype\\tn_barcodes\\tn_pairs" > {output}
        bcftools query -f '{wildcards.population}\\t%CHROM\\t%POS\\t%END\\t%SVLEN\\t%SVTYPE\\t%BARCODES\\t%PAIRS\\n' {input} >> {output}
        """

rule aggregate_variants:
    input:
        collect("reports/data/{population}.sv.stats", population = populations)
    output:
        "inversions.bedpe",
        "deletions.bedpe",
        "duplications.bedpe",
        "breakends.bedpe"
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

rule configure_report:
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

rule group_reports:
    input: 
        "reports/_quarto.yml",
        "reports/_harpy.scss",
        faidx     = f"{workflow_geno}.fai",
        statsfile = "reports/data/{population}.sv.stats",
        qmd       = "workflow/report/leviathan.qmd"
    output:
        report = "reports/{population}.leviathan.html",
        qmd = temp("reports/{population}.leviathan.qmd")
    log:
        "logs/reports/{population}.report.log"
    params:
        sample= lambda wc: "-P sample:" + wc.get('population'),
        contigs= f"-P contigs:{plot_contigs}"
    conda:
        "envs/r.yaml"
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        FAIDX=$(realpath {input.faidx})
        STATS=$(realpath {input.statsfile})
        quarto render {output.qmd} --log {log} --quiet -P faidx:$FAIDX -P statsfile:$STATS {params}
        """

rule aggregate_report:
    input: 
        "reports/_quarto.yml",
        "reports/_harpy.scss",
        faidx      = f"{workflow_geno}.fai",
        statsfiles = collect("reports/data/{pop}.sv.stats", pop = populations),
        qmd        = "workflow/report/leviathan_pop.qmd"
    output:
        report = "reports/leviathan.summary.html",
        qmd = temp("reports/leviathan.summary.qmd")
    log:
        "logs/reports/summary.report.log"
    params:
        statsdir = "reports/data/",
        contigs = f"-P contigs:{plot_contigs}"
    conda:
        "envs/r.yaml"
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        FAIDX=$(realpath {input.faidx})
        INPATH=$(realpath {params.statsdir})
        quarto render {output.qmd} --log {log} --quiet -P faidx:$FAIDX -P statsdir:$INPATH {params.contigs}
        """

rule workflow_summary:
    default_target: True
    input:
        vcf = collect("vcf/{pop}.bcf", pop = populations),
        bedpe_agg = collect("{sv}.bedpe", sv = ["inversions", "deletions","duplications", "breakends"]),
        reports = collect("reports/{pop}.leviathan.html", pop = populations) if not skip_reports else [],
        agg_report = "reports/leviathan.summary.html" if not skip_reports else []
    params:
        min_size = f"-v {min_size}",
        min_bc = f"-c {min_bc}",
        iters  = f"-B {iterations}",
        extra = extra
    run:
        summary = ["The harpy sv leviathan workflow ran using these parameters:"]
        summary.append(f"The provided reference genome: {bn}")
        concat = "The alignments were concatenated using:\n"
        concat += "\tconcatenate_bam.py --bx -o groupname.bam -b samples.list"
        summary.append(concat)
        bc_idx = "The barcodes were indexed using:\n"
        bc_idx += "LRez index bam -p -b INPUT"
        summary.append(bc_idx)
        svcall = "Leviathan was called using:\n"
        svcall += f"\tLEVIATHAN -b INPUT -i INPUT.BCI -g GENOME {params}"
        summary.append(svcall)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['snakemake']['relative']}"
        summary.append(sm)
        with open("workflow/sv.leviathan.summary", "w") as f:
            f.write("\n\n".join(summary))