containerized: "docker://pdimens/harpy:latest"

import os
import sys
import subprocess
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

##TODO MANUAL PRUNING OF SWITCH ERRORS
# https://github.com/vibansal/HapCUT2/blob/master/outputformat.md

pruning           = config["prune"]
molecule_distance = config["molecule_distance"]
extra             = config.get("extra", "") 
outdir 			  = config["output_directory"]
envdir            = os.getcwd() + "/.harpy_envs"
skipreports       = config["skip_reports"]
samples_from_vcf  = config["samples_from_vcf"]
variantfile       = config["inputs"]["variantfile"]
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
    geno       = f"Genome/{bn}.fasta"
    genofai    = f"Genome/{bn}.fasta.fai"
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
        outdir + "/workflow/input/vcf/{sample}.het.vcf"
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
        outdir + "/workflow/input/vcf/{sample}.hom.vcf"
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
        vcf = outdir + "/workflow/input/vcf/{sample}.het.vcf",
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
    benchmark:
        ".Benchmark/Phase/extracthairs.{sample}.txt"
    shell:
        "extractHAIRS {params} --nf 1 --bam {input.bam} --VCF {input.vcf} --out {output} > {log} 2>&1"

rule link_fragments:
    input: 
        bam       = get_alignments,
        vcf       = outdir + "/workflow/input/vcf/{sample}.het.vcf",
        fragments = outdir + "/extract_hairs/{sample}.unlinked.frags"
    output:
        outdir + "/link_fragments/{sample}.linked.frags"
    log:
        outdir + "/logs/link_fragments/{sample}.linked.log"
    params:
        d = molecule_distance
    conda:
        f"{envdir}/phase.yaml"
    benchmark:
        ".Benchmark/Phase/linkfrag.{sample}.txt"
    shell:
        "LinkFragments.py --bam {input.bam} --VCF {input.vcf} --fragments {input.fragments} --out {output} -d {params} > {log} 2>&1"

rule phase:
    input:
        vcf       = outdir + "/workflow/input/vcf/{sample}.het.vcf",
        fragments = fragfile
    output: 
        blocks    = outdir + "/phase_blocks/{sample}.blocks",
        vcf       = temp(outdir + "/phase_blocks/{sample}.blocks.phased.VCF")
    log:
        outdir + "/logs/hapcut2/{sample}.blocks.phased.log"
    params: 
        prune = f"--threshold {pruning}" if pruning > 0 else "--no_prune 1",
        extra = extra
    conda:
        f"{envdir}/phase.yaml"
    benchmark:
        ".Benchmark/Phase/phase.{sample}.txt"
    shell:
        "HAPCUT2 --fragments {input.fragments} --vcf {input.vcf} {params} --out {output.blocks} --nf 1 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1 > {log} 2>&1"

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
        outdir + "/workflow/input/vcf/{sample}.hom.vcf"
    output:
        outdir + "/workflow/input/gzvcf/{sample}.hom.vcf.gz"

rule merge_het_hom:
    input:
        phase = outdir + "/phase_blocks/{sample}.phased.vcf.gz",
        orig  = outdir + "/workflow/input/gzvcf/{sample}.hom.vcf.gz"
    output:
        bcf = outdir + "/phased_samples/{sample}.phased.annot.bcf",
        idx = outdir + "/phased_samples/{sample}.phased.annot.bcf.csi"
    threads:
        2
    benchmark:
        ".Benchmark/Phase/mergeAnno.{sample}.txt"
    container:
        None
    shell:
        "bcftools annotate -Ob -o {output.bcf} --write-index -a {input.phase} -c CHROM,POS,FMT/GT,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT {input.orig}"

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
            parse_phaseblocks.py -i $i >> {params}
        done
        gzip {params}
        """

rule phase_report:
    input:
        outdir + "/reports/blocks.summary.gz"
    output:
        outdir + "/reports/phase.html"
    log:
        logfile = outdir + "/logs/report.log"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/hapcut.Rmd"

rule workflow_summary:
    default_target: True
    input:
        vcf = outdir + "/variants.phased.bcf",
        reports = outdir + "/reports/phase.html" if not skipreports else []
    params:
        prune = f"--threshold {pruning}" if pruning > 0 else "--no_prune 1",
        extra = extra
    run:
        import glob
        for logfile in glob.glob(f"{outdir}/logs/**/*", recursive = True):
            if os.path.isfile(logfile) and os.path.getsize(logfile) == 0:
                os.remove(logfile)
        for logfile in glob.glob(f"{outdir}/logs/**/*", recursive = True):
            if os.path.isdir(logfile) and not os.listdir(logfile):
                os.rmdir(logfile)
        with open(outdir + "/workflow/phase.summary", "w") as f:
            _ = f.write("The harpy phase workflow ran using these parameters:\n\n")
            _ = f.write(f"The provided variant file: {variantfile}\n")
            _ = f.write("The variant file was split by sample and filtered for heterozygous sites using:\n")
            _ = f.write("""    bcftools view -s SAMPLE | bcftools view -i 'GT="het"' \n""")
            _ = f.write("Phasing was performed using the components of HapCut2:\n")
            _ = f.write("    extractHAIRS " + linkarg + " --nf 1 --bam sample.bam --VCF sample.vcf --out sample.unlinked.frags\n")
            _ = f.write("    LinkFragments.py --bam sample.bam --VCF sample.vcf --fragments sample.unlinked.frags --out sample.linked.frags -d " + f"{molecule_distance}" + "\n")
            _ = f.write("    HAPCUT2 --fragments sample.linked.frags --vcf sample.vcf --out sample.blocks --nf 1 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1" + f" {params[0]} {params[1]}" + "\n\n")
            _ = f.write("Variant annotation was performed using:\n")
            _ = f.write("    bcftools annotate -a sample.phased.vcf -c CHROM,POS,FMT/GT,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT \n")
            _ = f.write("    bcftools merge --output-type b samples.annot.bcf\n\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")
