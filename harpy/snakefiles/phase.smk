containerized: "docker://pdimens/harpy:latest"

import os
import subprocess
import logging
from pathlib import Path

onstart:
    logger.addHandler(logging.FileHandler(config["snakemake_log"]))
onsuccess:
    os.remove(logger.logfile)
onerror:
    os.remove(logger.logfile)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

outdir 			  = config["output_directory"]
workflowdir       = f"{outdir}/workflow"
pruning           = config["prune"]
molecule_distance = config["molecule_distance"]
extra             = config.get("extra", "") 
envdir            = os.path.join(os.getcwd(), outdir, "workflow", "envs")
samples_from_vcf  = config["samples_from_vcf"]
variantfile       = config["inputs"]["variantfile"]
skip_reports      = config["reports"]["skip"]
plot_contigs      = config["reports"]["plot_contigs"]
bamlist     = config["inputs"]["alignments"]
bamdict     = dict(zip(bamlist, bamlist))
if config["ignore_bx"]:
    fragfile = outdir + "/extract_hairs/{sample}.unlinked.frags"
    linkarg = "--10x 0"
else:
    fragfile =  outdir + "/link_fragments/{sample}.linked.frags"
    linkarg  = "--10x 1"
if samples_from_vcf:
    bcfquery = subprocess.Popen(["bcftools", "query", "-l", variantfile], stdout=subprocess.PIPE)
    samplenames = bcfquery.stdout.read().decode().split()
else:
    samplenames = [Path(i).stem for i in bamlist]
if config["inputs"].get("genome", None):
    genomefile = config["inputs"]["genome"]
    if genomefile.lower().endswith(".gz"):
        bn = Path(Path(genomefile).stem).stem
    else:
        bn = Path(genomefile).stem
    geno       = f"{workflowdir}/genome/{bn}"
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
        workflowdir + "/input/vcf/{sample}.het.vcf"
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
        workflowdir + "/input/vcf/{sample}.hom.vcf"
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
    rule process_genome:
        input:
            genomefile
        output: 
            geno
        container:
            None
        shell: 
            "seqtk seq {input} > {output}"

    rule index_genome:
        input: 
            geno
        output: 
            genofai
        log:
            f"Genome/{bn}.faidx.log"
        container:
            None
        shell: 
            "samtools faidx --fai-idx {output} {input} 2> {log}"

rule extract_hairs:
    input:
        vcf = workflowdir + "/input/vcf/{sample}.het.vcf",
        bam = get_alignments,
        bai = get_align_index,
        geno = geno,
        fai  = genofai
    output:
        outdir + "/extract_hairs/{sample}.unlinked.frags"
    log:
        outdir + "/logs/extract_hairs/{sample}.unlinked.log"
    params:
        indels = indelarg,
        bx = linkarg
    conda:
        f"{envdir}/phase.yaml"
    shell:
        "extractHAIRS {params} --nf 1 --bam {input.bam} --VCF {input.vcf} --out {output} > {log} 2>&1"

rule link_fragments:
    input: 
        bam       = get_alignments,
        vcf       = workflowdir + "/input/vcf/{sample}.het.vcf",
        fragments = outdir + "/extract_hairs/{sample}.unlinked.frags"
    output:
        outdir + "/link_fragments/{sample}.linked.frags"
    log:
        outdir + "/logs/link_fragments/{sample}.linked.log"
    params:
        d = molecule_distance
    conda:
        f"{envdir}/phase.yaml"
    shell:
        "LinkFragments.py --bam {input.bam} --VCF {input.vcf} --fragments {input.fragments} --out {output} -d {params} > {log} 2>&1"

rule phase:
    input:
        vcf       = workflowdir + "/input/vcf/{sample}.het.vcf",
        fragments = fragfile
    output: 
        blocks    = outdir + "/phase_blocks/{sample}.blocks",
        vcf       = temp(outdir + "/phase_blocks/{sample}.blocks.phased.VCF")
    log:
        outdir + "/logs/hapcut2/{sample}.blocks.phased.log"
    params: 
        prune = f"--threshold {pruning}" if pruning > 0 else "--no_prune 1",
        fixed_params = "--nf 1 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1",
        extra = extra
    conda:
        f"{envdir}/phase.yaml"
    shell:
        "HAPCUT2 --fragments {input.fragments} --vcf {input.vcf} --out {output.blocks} {params} > {log} 2>&1"

rule compress_phaseblock:
    input:
        outdir + "/phase_blocks/{sample}.blocks.phased.VCF"
    output:
        outdir + "/phase_blocks/{sample}.phased.vcf.gz"
    container:
        None
    shell:
        "bcftools view -Oz6 -o {output} --write-index {input}"

use rule compress_phaseblock as compress_vcf with:
    input:
        workflowdir + "/input/vcf/{sample}.hom.vcf"
    output:
        workflowdir + "/input/gzvcf/{sample}.hom.vcf.gz"

rule merge_het_hom:
    priority: 100
    input:
        phase = outdir + "/phase_blocks/{sample}.phased.vcf.gz",
        orig  = workflowdir + "/input/gzvcf/{sample}.hom.vcf.gz"
    output:
        bcf = outdir + "/phased_samples/{sample}.phased.annot.bcf",
        idx = outdir + "/phased_samples/{sample}.phased.annot.bcf.csi"
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
        bcf = collect(outdir + "/phased_samples/{sample}.phased.annot.bcf", sample = samplenames),
        idx = collect(outdir + "/phased_samples/{sample}.phased.annot.bcf.csi", sample = samplenames)
    output:
        bcf = outdir + "/variants.phased.bcf",
        idx = outdir + "/variants.phased.bcf.csi"
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
        collect(outdir + "/phase_blocks/{sample}.blocks", sample = samplenames)
    output:
        outdir + "/reports/blocks.summary.gz"
    params:
        outdir + "/reports/blocks.summary"
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

rule report_config:
    input:
        yaml = f"{workflowdir}/report/_quarto.yml",
        scss = f"{workflowdir}/report/_harpy.scss"
    output:
        yaml = temp(f"{outdir}/reports/_quarto.yml"),
        scss = temp(f"{outdir}/reports/_harpy.scss")
    run:
        import shutil
        for i,o in zip(input,output):
            shutil.copy(i,o)

rule phase_report:
    input:
        f"{outdir}/reports/_quarto.yml",
        f"{outdir}/reports/_harpy.scss",
        data = f"{outdir}/reports/blocks.summary.gz",
        qmd = f"{workflowdir}/report/hapcut.qmd"
    output:
        html = f"{outdir}/reports/phase.html",
        qmd = temp(f"{outdir}/reports/phase.qmd")
    log:
        f"{outdir}/logs/report.log"
    params:
        f"-P contigs:{plot_contigs}"
    conda:
        f"{envdir}/r.yaml"
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        INFILE=$(realpath {input.data})
        quarto render {output.qmd} --log {log} --quiet -P blockfile:$INFILE
        """

rule workflow_summary:
    default_target: True
    input:
        vcf = outdir + "/variants.phased.bcf",
        reports = outdir + "/reports/phase.html" if not skip_reports else []
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
        sm = f"\t{config['workflow_call']}"
        summary.append(sm)
        with open(f"{workflowdir}/phase.summary", "w") as f:
            f.write("\n\n".join(summary))
