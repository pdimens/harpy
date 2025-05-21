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

pruning           = config["prune"]
molecule_distance = config["barcodes"]["distance_threshold"]
extra             = config.get("extra", "") 
samples_from_vcf  = config["samples_from_vcf"]
variantfile       = config["inputs"]["variantfile"]
skip_reports      = config["reports"]["skip"]
plot_contigs      = config["reports"]["plot_contigs"]
bamlist     = config["inputs"]["alignments"]
bamdict     = dict(zip(bamlist, bamlist))
if config["barcodes"]["ignore"]:
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

rule extract_het:
    input: 
        vcf = variantfile
    output:
        "worklfow/input/vcf/{sample}.het.vcf"
    container:
        None
    shell:
        """
        bcftools view -s {wildcards.sample} -Ou -m 2 -M 2 {input.vcf} | bcftools view -i 'GT="het"' > {output}
        """

rule extract_hom:
    input: 
        vcf = variantfile
    output:
        "worklfow/input/vcf/{sample}.hom.vcf"
    container:
        None
    shell:
        """
        bcftools view -s {wildcards.sample} -Ou {input.vcf} | bcftools view -i 'GT="hom"' > {output}
        """

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
            fasta = geno,
            fai = genofai
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
    input:
        vcf = "worklfow/input/vcf/{sample}.het.vcf",
        bam = get_alignments,
        bai = get_align_index,
        geno = geno,
        fai  = genofai
    output:
        "extract_hairs/{sample}.unlinked.frags"
    log:
        "logs/extract_hairs/{sample}.unlinked.log"
    params:
        indels = indelarg,
        bx = linkarg
    conda:
        "envs/phase.yaml"
    shell:
        "extractHAIRS {params} --nf 1 --bam {input.bam} --VCF {input.vcf} --out {output} > {log} 2>&1"

rule link_fragments:
    input: 
        bam       = get_alignments,
        vcf       = "worklfow/input/vcf/{sample}.het.vcf",
        fragments = "extract_hairs/{sample}.unlinked.frags"
    output:
        "link_fragments/{sample}.linked.frags"
    log:
        "logs/link_fragments/{sample}.linked.log"
    params:
        d = molecule_distance
    conda:
        "envs/phase.yaml"
    shell:
        "LinkFragments.py --bam {input.bam} --VCF {input.vcf} --fragments {input.fragments} --out {output} -d {params} > {log} 2>&1"

rule phase:
    input:
        vcf       = "worklfow/input/vcf/{sample}.het.vcf",
        fragments = fragfile
    output: 
        blocks    = "phase_blocks/{sample}.blocks",
        vcf       = temp("phase_blocks/{sample}.blocks.phased.VCF")
    log:
        "logs/hapcut2/{sample}.blocks.phased.log"
    params: 
        prune = f"--threshold {pruning}" if pruning > 0 else "--no_prune 1",
        fixed_params = "--nf 1 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1",
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

use rule compress_phaseblock as compress_vcf with:
    input:
        "worklfow/input/vcf/{sample}.hom.vcf"
    output:
        "worklfow/input/gzvcf/{sample}.hom.vcf.gz"

rule merge_het_hom:
    priority: 100
    input:
        phase = "phase_blocks/{sample}.phased.vcf.gz",
        orig  = "worklfow/input/gzvcf/{sample}.hom.vcf.gz"
    output:
        bcf = "phased_samples/{sample}.phased.annot.bcf",
        idx = "phased_samples/{sample}.phased.annot.bcf.csi"
    params:
        "-Ob --write-index -c CHROM,POS,FMT/GT,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT"
    threads:
        2
    container:
        None
    shell:
        "bcftools annotate -a {input.phase} -o {output.bcf} {params} {input.orig}"

rule merge_samples:
    input: 
        bcf = collect("phased_samples/{sample}.phased.annot.bcf", sample = samplenames),
        idx = collect("phased_samples/{sample}.phased.annot.bcf.csi", sample = samplenames)
    output:
        bcf = "variants.phased.bcf",
        idx = "variants.phased.bcf.csi"
    params:
        "true" if len(samplenames) > 1 else "false"
    threads:
        workflow.cores
    container:
        None
    shell:
        """
        if [ "{params}" = true ]; then
            bcftools merge --threads {threads} -Ob -o {output.bcf} --write-index {input.bcf}
        else
           cp {input.bcf} {output.bcf}
           cp {input.idx} {output.idx}
        fi
        """

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
            parse_phaseblocks.py $i >> {params}
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
        "envs/r.yaml"
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        INFILE=$(realpath {input.data})
        quarto render {output.qmd} --log {log} --quiet -P blockfile:$INFILE
        """

rule workflow_summary:
    default_target: True
    input:
        vcf = "variants.phased.bcf",
        reports = "reports/phase.html" if not skip_reports else []
    params:
        prune = f"--threshold {pruning}" if pruning > 0 else "--no_prune 1",
        extra = extra
    run:
        summary = ["The harpy phase workflow ran using these parameters:"]
        summary.append(f"The provided variant file: {variantfile}")
        hetsplit = "The variant file was split by sample and filtered for heterozygous sites using:\n"
        hetsplit += "\tbcftools view -s SAMPLE | bcftools view -i \'GT=\"het\"\'"
        summary.append(hetsplit)
        phase = "Phasing was performed using the components of HapCut2:\n"
        phase += "\textractHAIRS {linkarg} --nf 1 --bam sample.bam --VCF sample.vcf --out sample.unlinked.frags\n"
        phase += f"\tLinkFragments.py --bam sample.bam --VCF sample.vcf --fragments sample.unlinked.frags --out sample.linked.frags -d {molecule_distance}\n"
        phase += f"\tHAPCUT2 --fragments sample.linked.frags --vcf sample.vcf --out sample.blocks --nf 1 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1 {params.prune} {params.extra}\n"
        summary.append(phase)
        annot = "Variant annotation was performed using:\n"
        annot += "\tbcftools annotate -a sample.phased.vcf -c CHROM,POS,FMT/GT,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT\n"
        annot += "\tbcftools merge --output-type b samples.annot.bcf"
        summary.append(annot)
        sm = "The Snakemake workflow was called via command line:\n"
        sm = f"\t{config['snakemake']['relative']}"
        summary.append(sm)
        with open("workflow/phase.summary", "w") as f:
            f.write("\n\n".join(summary))
