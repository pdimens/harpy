containerized: "docker://pdimens/harpy:latest"

import os
import subprocess
import logging
from pathlib import Path

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake"]["log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

bc_type           = config["linkedreads"]["type"]
pruning           = config["phasing"]["prune"]
map_qual          = config["phasing"]["min_map_quality"]
base_qual         = config["phasing"]["min_base_quality"]
molecule_distance = config["linkedreads"]["distance_threshold"]
extra             = config.get("extra", "") 
samples_from_vcf  = config["inputs"]["vcf"]["prioritize_samples"]
variantfile       = config["inputs"]["vcf"]["file"]
skip_reports      = config["reports"]["skip"]
plot_contigs      = config["reports"]["plot_contigs"]
bamlist           = config["inputs"]["alignments"]
bamdict           = dict(zip(bamlist, bamlist))
invalid_regex = {
    "haplotagging" : "'$4 !~ /[ABCD]00/'",
    "stlfr" : "'$4 !~ /^0_|_0_|_0$/'",
    "tellseq": "'$4 !~ /N/'"
}
if bc_type == "none":
    fragfile = "extract_hairs/{sample}.unlinked.frags"
    linkarg = "--10x 0"
else:
    fragfile =  "link_fragments/{sample}.linked.frags"
    linkarg  = "--10x 1"
if samples_from_vcf:
    bcfquery = subprocess.Popen(["bcftools", "query", "-l", variantfile], stdout=subprocess.PIPE)
    samplenames = bcfquery.stdout.read().decode().split()
else:
    samplenames = [Path(i).stem for i in bamlist]
if config["inputs"].get("reference", None):
    genomefile = config["inputs"]["reference"]
    if genomefile.lower().endswith(".gz"):
        bn = Path(Path(genomefile).stem).stem
    else:
        bn = Path(genomefile).stem
    geno       = f"workflow/reference/{bn}"
    genofai    = f"{geno}.fai"
    indelarg   = f"--indels 1 --ref {geno}"
    indels     = True
else:
    indelarg   = ""
    genomefile = []
    geno       = []
    genofai    = []
    bn         = []
    indels     = False

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

rule isolate_sample:
    input: 
        variantfile
    output:
        vcf = temp("workflow/input/original/{sample}.bcf"),
        csi = temp("workflow/input/original/{sample}.bcf.csi")
    container:
        None
    shell:
        "bcftools view -Ob -W -s {wildcards.sample} -o {output.vcf} {input}"

rule isolate_het_snps:
    input: 
        "workflow/input/original/{sample}.bcf"
    output:
        temp("workflow/input/heterozygotes/{sample}.het.vcf")
    container:
        None
    shell:
        "bcftools view -m 2 -M 2 -i 'GT=\"het\"' {input} > {output}"

rule index_alignments:
    input:
        lambda wc: bamdict[wc.bam]
    output:
        "{bam}.bai"
    container:
        None
    shell:
        "samtools index {input}"

if indels:
    rule preprocess_reference:
        input:
            genomefile
        output: 
            fasta = temp(geno),
            fai = temp(genofai)
        log:
            f"workflow/reference/{bn}.preprocess.log"
        container:
            None
        shell: 
            """
            seqtk seq {input} > {output.fasta}
            samtools faidx --fai-idx {output.fai} {output.fasta} 2> {log}
            """

rule extract_hairs:
    retries:
        2
    input:
        vcf = "workflow/input/heterozygotes/{sample}.het.vcf",
        bam = get_alignments,
        bai = get_align_index,
        geno = geno,
        fai  = genofai
    output:
        all_bc = "extract_hairs/{sample}.unlinked.all.frags",
        no_invalid = "extract_hairs/{sample}.unlinked.frags"
    log:
        "logs/extract_hairs/{sample}.unlinked.log"
    params:
        static = f"{indelarg} {linkarg} --mmq {map_qual} --mbq {base_qual} --nf 1 --maxfragments 1500000",
        purge_invalid = invalid_regex.get(bc_type, "'$4 !~ /N/'")
    conda:
        "envs/phase.yaml"
    shell:
        """
        extractHAIRS {params.static} --bam {input.bam} --VCF {input.vcf} --out {output.all_bc} > {log} 2>&1
        awk {params.purge_invalid} {output.all_bc} > {output.no_invalid}
        """

rule link_fragments:
    input: 
        bam       = get_alignments,
        vcf       = "workflow/input/heterozygotes/{sample}.het.vcf",
        fragments = "extract_hairs/{sample}.unlinked.frags"
    output:
        "link_fragments/{sample}.linked.frags"
    log:
        "logs/link_fragments/{sample}.linked.log"
    params:
        f"-d {molecule_distance} --use-tag"
    conda:
        "envs/phase.yaml"
    shell:
        "LinkFragments.py --bam {input.bam} --VCF {input.vcf} --fragments {input.fragments} --out {output} {params} > {log} 2>&1"

rule phase:
    input:
        vcf       = "workflow/input/heterozygotes/{sample}.het.vcf",
        fragments = fragfile
    output: 
        blocks    = "phase_blocks/{sample}.blocks",
        vcf       = temp("phase_blocks/{sample}.blocks.phased.VCF")
    log:
        "logs/hapcut2/{sample}.blocks.phased.log"
    params: 
        prune = f"--threshold {pruning}" if pruning > 0 else "--no_prune 1",
        fixed_params = "--nf 1 --error_analysis_mode 1 --outvcf 1 --call_homozygous 1",
        extra = extra
    conda:
        "envs/phase.yaml"
    shell:
        "HAPCUT2 --fragments {input.fragments} --vcf {input.vcf} --out {output.blocks} {params} > {log} 2>&1"

rule compress_phaseblock:
    input:
        "phase_blocks/{sample}.blocks.phased.VCF"
    output:
        "phase_blocks/{sample}.phased.vcf.gz"
    container:
        None
    shell:
        "bcftools view -Oz6 -o {output} --write-index {input}"

rule annotate_phase:
    priority: 100
    input:
        "workflow/input/original/{sample}.bcf.csi",
        phase = "phase_blocks/{sample}.phased.vcf.gz",
        orig = "workflow/input/original/{sample}.bcf"
    output:
        bcf = "phased_samples/{sample}.phased.annot.bcf",
        idx = "phased_samples/{sample}.phased.annot.bcf.csi"
    log:
        "logs/annotate/{sample}.annotate.log"
    params:
        "-Ob --write-index -c CHROM,POS,FMT/GT,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT"
    threads:
        2
    container:
        None
    shell:
        "bcftools annotate -a {input.phase} -o {output.bcf} {params} {input.orig} 2> {log}"
rule create_merge_list:
    input:
        bcf = collect("phased_samples/{sample}.phased.annot.bcf", sample = samplenames)
    output:
        filelist = temp("phased_samples/bcf.list")
    run:
        with open(output.filelist, "w") as f_out:
            for i in input.bcf:
                f_out.write(i + "\n")

rule merge_samples:
    priority: 100
    input: 
        collect("phased_samples/{sample}.phased.annot.{ext}", sample = samplenames, ext = ["bcf", "bcf.csi"]),
        filelist = "phased_samples/bcf.list"
    output:
        "variants.phased.bcf.csi",
        bcf = "variants.phased.bcf"
    threads:
        workflow.cores
    container:
        None
    shell:
        "bcftools merge --threads {threads} --force-single -l {input.filelist} -Ob -o {output.bcf} --write-index"

rule summarize_blocks:
    input:
        collect("phase_blocks/{sample}.blocks", sample = samplenames)
    output:
        "reports/blocks.summary.gz"
    params:
        "reports/blocks.summary"
    container:
        None
    shell:
        """
        echo -e "sample\\tcontig\\tn_snp\\tpos_start\\tblock_length" > {params}
        for i in {input}; do
            parse_phaseblocks $i >> {params}
        done
        gzip {params}
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

rule phase_report:
    input:
        "reports/_quarto.yml",
        "reports/_harpy.scss",
        data = "reports/blocks.summary.gz",
        qmd = "workflow/report/hapcut.qmd"
    output:
        html = "reports/phase.html",
        qmd = temp("reports/phase.qmd")
    log:
        "logs/report.log"
    params:
        f"-P contigs:{plot_contigs}"
    conda:
        "envs/report.yaml"
    retries:
        3
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        INFILE=$(realpath {input.data})
        quarto render {output.qmd} --no-cache --log {log} --quiet -P blockfile:$INFILE
        """

rule workflow_summary:
    default_target: True
    input:
        vcf = "variants.phased.bcf",
        reports = "reports/phase.html" if not skip_reports else []
