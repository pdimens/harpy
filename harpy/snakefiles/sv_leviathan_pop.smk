containerized: "docker://pdimens/harpy:latest"

import os
import re
import logging

onstart:
    logger.logger.addHandler(logging.FileHandler(config["snakemake_log"]))
onsuccess:
    os.remove(logger.logfile)
onerror:
    os.remove(logger.logfile)
wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+",
    population = "[a-zA-Z0-9._-]+"

envdir      = os.getcwd() + "/.harpy_envs"
genomefile 	= config["inputs"]["genome"]
bamlist     = config["inputs"]["alignments"]
groupfile 	= config["inputs"]["groupings"]
extra 		= config.get("extra", "") 
min_sv      = config["min_sv"]
min_bc      = config["min_barcodes"]
iterations  = config["iterations"]
outdir      = config["output_directory"]
skip_reports = config["skip_reports"]
bn 			= os.path.basename(genomefile)
if bn.lower().endswith(".gz"):
    bn = bn[:-3]
    
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

rule preproc_groups:
    input:
        groupfile
    output:
        outdir + "/workflow/sample.groups"
    run:
        with open(input[0], "r") as infile, open(output[0], "w") as outfile:
            _ = [outfile.write(i) for i in infile.readlines() if not i.lstrip().startswith("#")]

rule concat_list:
    input:
        outdir + "/workflow/sample.groups"
    output:
        outdir + "/workflow/merge_samples/{population}.list"
    run:
        with open(output[0], "w") as fout:
            for bamfile in popdict[wildcards.population]:
                _ = fout.write(bamfile + "\n")

rule concat_groups:
    input: 
        bamlist  = outdir + "/workflow/merge_samples/{population}.list",
        bamfiles = lambda wc: collect("{sample}", sample = popdict[wc.population]) 
    output:
        temp(outdir + "/workflow/input/{population}.unsort.bam")
    log:
        outdir + "/logs/concat_groups/{population}.concat.log"
    threads:
        1
    container:
        None
    shell:
        "concatenate_bam.py -o {output} -b {input.bamlist} 2> {log}"

rule sort_groups:
    input:
        outdir + "/workflow/input/{population}.unsort.bam"
    output:
        bam = temp(outdir + "/workflow/input/{population}.bam"),
        bai = temp(outdir + "/workflow/input/{population}.bam.bai")
    log:
        outdir + "/logs/samtools_sort/{population}.sort.log"
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
        bam = outdir + "/workflow/input/{population}.bam",
        bai = outdir + "/workflow/input/{population}.bam.bai"
    output:
        temp(outdir + "/lrez_index/{population}.bci")
    benchmark:
        ".Benchmark/leviathan-pop/{population}.lrez"
    threads:
        max(10, workflow.cores)
    conda:
        f"{envdir}/variants.yaml"
    shell:
        "LRez index bam -p -b {input.bam} -o {output} --threads {threads}"

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
        bam    = outdir + "/workflow/input/{population}.bam",
        bai    = outdir + "/workflow/input/{population}.bam.bai",
        bc_idx = outdir + "/lrez_index/{population}.bci",
        genome = f"Genome/{bn}",
        genidx = multiext(f"Genome/{bn}", ".fai", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output:
        vcf = temp(outdir + "/vcf/{population}.vcf"),
        candidates = outdir + "/logs/leviathan/{population}.candidates"
    log:  
        runlog = outdir + "/logs/leviathan/{population}.leviathan.log",
    params:
        min_sv = f"-v {min_sv}",
        min_bc = f"-c {min_bc}",
        iters  = f"-B {iterations}",
        extra = extra
    threads:
        workflow.cores - 1
    conda:
        f"{envdir}/variants.yaml"
    benchmark:
        ".Benchmark/leviathan-pop/{population}.variantcall"
    shell:
        "LEVIATHAN -b {input.bam} -i {input.bc_idx} {params} -g {input.genome} -o {output.vcf} -t {threads} --candidates {output.candidates} 2> {log.runlog}"

rule sort_variants:
    input:
        outdir + "/vcf/{population}.vcf"
    output:
        outdir + "/vcf/{population}.bcf"
    params:
        "{wildcards.population}"
    container:
        None
    shell:        
        "bcftools sort -Ob --output {output} {input} 2> /dev/null"

rule variant_stats:
    input: 
        outdir + "/vcf/{population}.bcf"
    output:
        temp(outdir + "/reports/data/{population}.sv.stats")
    container:
        None
    shell:
        """
        echo -e "population\\tcontig\\tposition_start\\tposition_end\\tlength\\ttype\\tn_barcodes\\tn_pairs" > {output}
        bcftools query -f '{wildcards.population}\\t%CHROM\\t%POS\\t%END\\t%SVLEN\\t%SVTYPE\\t%BARCODES\\t%PAIRS\\n' {input} >> {output}
        """

rule aggregate_variants:
    input:
        collect(outdir + "/reports/data/{population}.sv.stats", population = populations)
    output:
        outdir + "/inversions.bedpe",
        outdir + "/deletions.bedpe",
        outdir + "/duplications.bedpe",
        outdir + "/breakends.bedpe"
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

rule group_reports:
    input:	
        statsfile = outdir + "/reports/data/{population}.sv.stats",
        bcf       = outdir + "/vcf/{population}.bcf",
        faidx     = f"Genome/{bn}.fai"
    output:
        outdir + "/reports/{population}.sv.html"
    log:
        logfile = outdir + "/logs/reports/{population}.report.log"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/leviathan.Rmd"

rule aggregate_report:
    input:	
        faidx      = f"Genome/{bn}.fai",
        statsfiles = collect(outdir + "/reports/data/{pop}.sv.stats", pop = populations)
    output:
        outdir + "/reports/leviathan.summary.html"
    log:
        logfile = outdir + "/logs/reports/summary.report.log"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/leviathan_pop.Rmd"

rule workflow_summary:
    default_target: True
    input:
        vcf = collect(outdir + "/vcf/{pop}.bcf", pop = populations),
        bedpe_agg = collect(outdir + "/{sv}.bedpe", sv = ["inversions", "deletions","duplications", "breakends"]),
        reports = collect(outdir + "/reports/{pop}.sv.html", pop = populations) if not skip_reports else [],
        agg_report = outdir + "/reports/leviathan.summary.html" if not skip_reports else []
    params:
        min_sv = f"-v {min_sv}",
        min_bc = f"-c {min_bc}",
        iters  = f"-B {iterations}",
        extra = extra
    run:
        summary = ["The harpy sv leviathan workflow ran using these parameters:"]
        summary.append(f"The provided genome: {bn}")
        concat = "The alignments were concatenated using:\n"
        concat += "\tconcatenate_bam.py -o groupname.bam -b samples.list"
        summary.append(concat)
        bc_idx = "The barcodes were indexed using:\n"
        bc_idx += "LRez index bam -p -b INPUT"
        summary.append(bc_idx)
        svcall = "Leviathan was called using:\n"
        svcall += f"\tLEVIATHAN -b INPUT -i INPUT.BCI -g GENOME {params}"
        summary.append(svcall)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config["workflow_call"]}"
        summary.append(sm)
        with open(outdir + "/workflow/sv.leviathan.summary", "w") as f:
            f.write("\n\n".join(summary))