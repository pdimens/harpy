containerized: "docker://pdimens/harpy:latest"

import os
import re
import logging
from pathlib import Path

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake"]["log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

genomefile  = config["inputs"]["reference"]
bamlist     = config["inputs"]["alignments"]
bamdict     = dict(zip(bamlist, bamlist))
samplenames = {Path(i).stem for i in bamlist}
extra       = config.get("extra", None) 
mol_dist    = config["molecule_distance"]
min_size      = config["min_size"]
min_barcodes = config["min_barcodes"]
min_quality  = config["min_quality"]
bn          = os.path.basename(genomefile)
genome_zip  = True if bn.lower().endswith(".gz") else False
workflow_geno = f"workflow/reference/{bn}"
workflow_geno_idx = f"{workflow_geno}.gzi" if genome_zip else f"{workflow_geno}.fai"
skip_reports = config["reports"]["skip"]
plot_contigs = config["reports"]["plot_contigs"]    

def process_args(args):
    argsDict = {
        "min_mapq" : min_quality,
        "d"        : mol_dist,
        "min_sv"   : min_size,
        "k"        : min_barcodes
    }
    if args:
        words = [i for i in re.split(r"\s|=", args) if len(i) > 0]
        for i in zip(words[::2], words[1::2]):
            if "blacklist" in i or "candidates" in i:
                argsDict[i[0].lstrip("-")] = i[1]
    return argsDict

argdict = process_args(extra)

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

rule naibr_config:
    input:
        bam = get_alignments
    output:
        cfg = "workflow/input/{sample}.naibr"
    params:
        smp = lambda wc: wc.get("sample"),
        thd = min(10, workflow.cores)
    run:
        with open(output[0], "w") as conf:
            _ = conf.write(f"bam_file={input.bam}\n")
            _ = conf.write(f"prefix={params.smp}\n")
            _ = conf.write(f"outdir={params.smp}\n")
            _ = conf.write(f"threads={params.thd}\n")
            for i in argdict:
                _ = conf.write(f"{i}={argdict[i]}\n")

rule call_variants:
    input:
        bam   = get_alignments,
        bai   = get_align_index,
        conf  = "workflow/input/{sample}.naibr"
    output:
        bedpe = temp("{sample}/{sample}.bedpe"),
        refmt = temp("{sample}/{sample}.reformat.bedpe"),
        vcf   = temp("{sample}/{sample}.vcf"),
        log   = temp("{sample}/{sample}.log")
    log:
        "logs/naibr/{sample}.naibr.log"
    threads:
        min(10, workflow.cores -1)
    conda:
        "envs/variants.yaml"     
    shell:
        "naibr {input.conf} > {log} 2>&1 && rm -rf naibrlog"

rule infer_variants:
    input:
        bedpe = "{sample}/{sample}.bedpe",
        refmt = "{sample}/{sample}.reformat.bedpe",
        vcf   = "{sample}/{sample}.vcf"
    output:
        bedpe = "bedpe/{sample}.bedpe",
        refmt = "IGV/{sample}.reformat.bedpe",
        fail  = "bedpe/qc_fail/{sample}.fail.bedpe",
        vcf   = "vcf/{sample}.vcf" 
    container:
        None
    shell:
        """
        infer_sv.py {input.bedpe} -f {output.fail} > {output.bedpe}
        cp {input.refmt} {output.refmt}
        cp {input.vcf} {output.vcf}
        """

rule aggregate_variants:
    input:
        collect("bedpe/{sample}.bedpe", sample = samplenames)
    output:
        "inversions.bedpe",
        "deletions.bedpe",
        "duplications.bedpe"
    run:
        from pathlib import Path
        with open(output[0], "w") as inversions, open(output[1], "w") as deletions, open(output[2], "w") as duplications:
            header = ["Sample","Chr1","Break1","Chr2","Break2","SplitMolecules","DiscordantReads","Orientation","Haplotype","Score","PassFilter","SV"]
            _ = inversions.write("\t".join(header) + "\n")
            _ = deletions.write("\t".join(header) + "\n")
            _ = duplications.write("\t".join(header) + "\n")
            for varfile in input:
                samplename = Path(varfile).stem
                with open(varfile, "r") as f_in:
                    # read the header to skip it
                    f_in.readline()
                    # read the rest of it
                    while True:
                        line = f_in.readline()
                        if not line:
                            break
                        record = line.rstrip().split("\t")
                        if record[-1] == "inversion":
                            _ = inversions.write(f"{samplename}\t{line}")
                        elif record[-1] == "deletion":
                            _ = deletions.write(f"{samplename}\t{line}")
                        elif record[-1] == "duplication":
                            _ = duplications.write(f"{samplename}\t{line}")

rule preprocess_reference:
    input:
        genomefile
    output: 
        geno = workflow_geno,
        fai = f"{workflow_geno}.fai",
        gzi = f"{workflow_geno}.gzi" if genome_zip else []
    log:
        f"{workflow_geno}.preprocess.log"
    params:
        f"--gzi-idx {workflow_geno}.gzi" if genome_zip else ""
    container:
        None
    shell: 
        """
        if (file {input} | grep -q compressed ) ;then
            # is regular gzipped, needs to be BGzipped
            seqtk seq {input} | bgzip -c > {output.geno} 2> {log}
        else
            # if BZgipped or isn't compressed, just copied
            cp -f {input} {output.geno} 2> {log}
        fi
        samtools faidx {params} --fai-idx {output.fai} {output.geno} 2>> {log}
        """

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

rule sample_reports:
    input: 
        "reports/_quarto.yml",
        "reports/_harpy.scss",
        faidx = f"{workflow_geno}.fai",
        bedpe = "bedpe/{sample}.bedpe",
        qmd   = "workflow/report/naibr.qmd"
    output:
        report = "reports/{sample}.naibr.html",
        qmd = temp("reports/{sample}.naibr.qmd")
    log:
        "logs/reports/{sample}.report.log"
    params:
        sample= lambda wc: "-P sample:" + wc.get('sample'),
        contigs= f"-P contigs:{plot_contigs}"
    conda:
        "envs/r.yaml"
    retries:
        3
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        FAIDX=$(realpath {input.faidx})
        BEDPE=$(realpath {input.bedpe})
        quarto render {output.qmd} --log {log} --quiet -P faidx:$FAIDX -P bedpe:$BEDPE {params}
        """

rule workflow_summary:
    default_target: True
    input:
        bedpe = collect("bedpe/{sample}.bedpe", sample = samplenames),
        bedpe_agg = collect("{sv}.bedpe", sv = ["inversions", "deletions","duplications"]),
        reports =  collect("reports/{sample}.naibr.html", sample = samplenames) if not skip_reports else []
    run:
        #os.system(f"rm -rf naibrlog")
        summary = ["The harpy sv naibr workflow ran using these parameters:"]
        summary.append(f"The provided reference genome: {bn}")
        naibr = "naibr variant calling ran using these configurations:\n"
        naibr += "\tbam_file=BAMFILE\n"
        naibr += "\tprefix=PREFIX\n"
        naibr += "\toutdir=Variants/naibr/PREFIX\n"
        naibr += "\n\t".join([f"{k}={v}" for k,v in argdict.items()])
        summary.append(naibr)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['snakemake']['relative']}"
        summary.append(sm)
        with open("workflow/sv.naibr.summary", "w") as f:
            f.write("\n\n".join(summary))
