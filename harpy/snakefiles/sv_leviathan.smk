containerized: "docker://pdimens/harpy:latest"

import os
import re
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
genomefile  = config["inputs"]["genome"]
bamlist     = config["inputs"]["alignments"]
bamdict     = dict(zip(bamlist, bamlist))
samplenames = {Path(i).stem for i in bamlist}
min_sv      = config["min_sv"]
min_bc      = config["min_barcodes"]
iterations  = config["iterations"]
extra       = config.get("extra", "") 
skip_reports = config["reports"]["skip"]
plot_contigs = config["reports"]["plot_contigs"]    
bn          = os.path.basename(genomefile)
genome_zip  = True if bn.lower().endswith(".gz") else False
if genome_zip:
    bn = bn[:-3]

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
        lambda wc: bamdict[wc.bam]
    output:
        "{bam}.bai"
    container:
        None
    shell:
        "samtools index {input}"

rule index_barcode:
    input: 
        bam = get_alignments,
        bai = get_align_index
    output:
        temp(outdir + "/lrezIndexed/{sample}.bci")
    threads:
        max(10, workflow.cores)
    container:
        None
    shell:
        "LRez index bam --threads {threads} -p -b {input.bam} -o {output}"

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
        f"{envdir}/variants.yaml"
    shell:
        "LEVIATHAN -b {input.bam} -i {input.bc_idx} {params} -g {input.genome} -o {output.vcf} -t {threads} --candidates {output.candidates} 2> {log.runlog}"

rule sort_variants:
    priority: 100
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
        with (
            open(output[0], "w") as inversions,
            open(output[1], "w") as deletions,
            open(output[2], "w") as duplications,
            open(output[3], "w") as breakends
        ):
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

rule report_config:
    input:
        f"{outdir}/workflow/report/_quarto.yml"
    output:
        f"{outdir}/reports/_quarto.yml"
    shell:
        "cp {input} {output}"

rule sample_reports:
    input: 
        faidx     = f"Genome/{bn}.fai",
        statsfile = outdir + "/reports/data/{sample}.sv.stats",
        qmd       = f"{outdir}/workflow/report/leviathan.qmd",
        yml       = f"{outdir}/reports/_quarto.yml"
    output:
        report = outdir + "/reports/{sample}.leviathan.html",
        qmd = temp(outdir + "/reports/{sample}.leviathan.qmd")
    log:
        outdir + "/logs/reports/{sample}.report.log"
    params:
        sample= lambda wc: "-P sample:" + wc.get('sample'),
        contigs= f"-P contigs:{plot_contigs}"
    conda:
        f"{envdir}/r.yaml"
    shell:
        """
        cp {input.qmd} {output.qmd}
        FAIDX=$(realpath {input.faidx})
        STATS=$(realpath {input.statsfile})
        quarto render {output.qmd} --log {log} --quiet -P faidx:$FAIDX -P statsfile:$STATS {params}
        """

rule workflow_summary:
    default_target: True
    input: 
        vcf = collect(outdir + "/vcf/{sample}.bcf", sample = samplenames),
        bedpe_agg = collect(outdir + "/{sv}.bedpe", sv = ["inversions", "deletions","duplications", "breakends"]),
        reports = collect(outdir + "/reports/{sample}.leviathan.html", sample = samplenames) if not skip_reports else []
    params:
        min_sv = f"-v {min_sv}",
        min_bc = f"-c {min_bc}",
        iters  = f"-B {iterations}",
        extra = extra
    run:
        summary = ["The harpy sv leviathan workflow ran using these parameters:"]
        summary.append(f"The provided genome: {bn}")
        bc_idx = "The barcodes were indexed using:\n"
        bc_idx += "LRez index bam -p -b INPUT"
        summary.append(bc_idx)
        svcall = "Leviathan was called using:\n"
        svcall += f"\tLEVIATHAN -b INPUT -i INPUT.BCI -g GENOME {params}"
        summary.append(svcall)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['workflow_call']}"
        summary.append(sm)
        with open(outdir + "/workflow/sv.leviathan.summary", "w") as f:
            f.write("\n\n".join(summary))