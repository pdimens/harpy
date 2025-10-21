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
#bamdict     = dict(zip(bamlist, bamlist))
vcffile     = config["inputs"]["vcf"]
samplenames = {Path(i).stem for i in bamlist}
extra       = config.get("extra", None) 
mol_dist    = config["molecule_distance"]
min_quality  = config["min_quality"]
min_size      = config["min_size"]
min_barcodes = config["min_barcodes"]
plot_contigs = config["reports"]["plot_contigs"]    
skip_reports = config["reports"]["skip"]
bn          = os.path.basename(genomefile)
if bn.lower().endswith(".gz"):
    workflow_geno = f"workflow/reference/{bn[:-3]}"
else:
    workflow_geno = f"workflow/reference/{bn}"
if vcffile.lower().endswith("bcf"):
    vcfindex = vcffile + ".csi"
else:
    vcfindex = vcffile + ".tbi"

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

rule preprocess_reference:
    input:
        genomefile
    output: 
        geno = workflow_geno,
        fai = f"{workflow_geno}.fai"
    log:
        f"{workflow_geno}.preprocess.log"
    shell: 
        """
        {{
            seqtk seq {input} > {output.geno}
            samtools faidx --fai-idx {output.fai} {output.geno}
        }} > {log}
        """

rule index_alignments:
    input:
        get_alignments
    output:
        "{sample}.bai"
    shell:
        "samtools index {input}"

rule index_snps:
    input:
        vcffile
    output:
        vcffile + ".csi"
    shell:
        "bcftools index {input}"

rule index_snps_gz:
    input:
        vcffile
    output:
        vcffile + ".tbi"
    shell:
        "tabix {input}"

rule phase_alignments:
    input:
        get_align_index,
        vcfindex,
        f"{workflow_geno}.fai",
        vcf = vcffile,
        aln = get_alignments,
        ref = workflow_geno
    output:
        bam = "phasedbam/{sample}.bam",
        log = "logs/whatshap-haplotag/{sample}.phase.log"
    params:
        mol_dist
    conda:
        "envs/phase.yaml"
    container:
        "docker://pdimens/harpy:phase_latest"
    threads:
        4
    shell:
        "whatshap haplotag --sample {wildcards.sample} --linked-read-distance-cutoff {params} --ignore-read-groups --tag-supplementary --output-threads={threads} -o {output.bam} --reference {input.ref} {input.vcf} {input.aln} 2> {output.log}"

rule log_phasing:
    input:
        collect("logs/whatshap-haplotag/{sample}.phase.log", sample = samplenames)
    output:
        "logs/whatshap-haplotag.log"
    shell:
        """
        echo -e "sample\\ttotal_alignments\\tphased_alignments" > {output}
        for i in {input}; do
            SAMP=$(basename $i .phaselog)
            echo -e "${{SAMP}}\\t$(grep "Total alignments" $i)\\t$(grep "could be tagged" $i)" |
                sed 's/ \\+ /\\t/g' | cut -f1,3,5 >> {output}
        done
        """

rule naibr_config:
    input:
        bam = "phasedbam/{sample}.bam"
    output:
        cfg = "workflow/input/{sample}.naibr"
    params:
        smp = lambda wc: wc.get("sample"),
        thd = min(10, workflow.cores - 1)
    run:
        with open(output.cfg, "w") as conf:
            _ = conf.write(f"bam_file={input.bam}\n")
            _ = conf.write(f"prefix={params.smp}\n")
            _ = conf.write(f"outdir={params.smp}\n")
            _ = conf.write(f"threads={params.thd}\n")
            for i in argdict:
                _ = conf.write(f"{i}={argdict[i]}\n")

rule index_phased:
    input:
        "phasedbam/{sample}.bam"
    output:
        "phasedbam/{sample}.bam.bai"
    shell:
        "samtools index {input} {output} 2> /dev/null"

rule call_variants:
    input:
        bam   = "phasedbam/{sample}.bam",
        bai   = "phasedbam/{sample}.bam.bai",
        conf  = "workflow/input/{sample}.naibr"
    output:
        bedpe = temp("{sample}/{sample}.bedpe"),
        refmt = temp("{sample}/{sample}.reformat.bedpe"),
        vcf   = temp("{sample}/{sample}.vcf"),
        log   = temp("{sample}/{sample}.log")
    log:
        "logs/naibr/{sample}.naibr.log"
    threads:
        10
    conda:
        "envs/variants.yaml"
    container:
        "docker://pdimens/harpy:variants_latest"
    shell:
        "naibr {input.conf} > {log} 2>&1 && rm -rf naibrlog"

rule infer_variants:
    priority: 100
    input:
        bedpe = "{sample}/{sample}.bedpe",
        refmt = "{sample}/{sample}.reformat.bedpe",
        vcf   = "{sample}/{sample}.vcf"
    output:
        bedpe = "bedpe/{sample}.bedpe",
        refmt = "IGV/{sample}.reformat.bedpe",
        fail  = "bedpe/qc_fail/{sample}.fail.bedpe",
        vcf   = "vcf/{sample}.vcf" 
    shell:
        """
        infer_sv {input.bedpe} -f {output.fail} > {output.bedpe}
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
        "envs/report.yaml"
    container:
        "docker://pdimens/harpy:report_latest"
    retries:
        3
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        FAIDX=$(realpath {input.faidx})
        BEDPE=$(realpath {input.bedpe})
        quarto render {output.qmd} --no-cache --log {log} --quiet -P faidx:$FAIDX -P bedpe:$BEDPE {params}
        """

rule all:
    default_target: True
    input:
        bedpe = collect("bedpe/{sample}.bedpe", sample = samplenames),
        bedpe_agg = collect("{sv}.bedpe", sv = ["inversions", "deletions","duplications"]),
        phaselog = "logs/whatshap-haplotag.log",
        reports =  collect("reports/{sample}.naibr.html", sample = samplenames) if not skip_reports else []
