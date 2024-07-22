containerized: "docker://pdimens/harpy:latest"

import sys
import subprocess
import multiprocessing
from pathlib import Path
from rich import print as rprint
from rich.panel import Panel

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

# toggle linked-read aware mode
if config["ignore_bx"]:
    fragfile = outdir + "/extractHairs/{sample}.unlinked.frags"
    linkarg = "--10x 0"
else:
    fragfile =  outdir + "/linkFragments/{sample}.linked.frags"
    linkarg  = "--10x 1"

if samples_from_vcf:
    bcfquery = subprocess.Popen(["bcftools", "query", "-l", variantfile], stdout=subprocess.PIPE)
    samplenames = bcfquery.stdout.read().decode().split()
else:
    samplenames = [Path(i).stem for i in bamlist]


# toggle indel mode
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

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

onerror:
    print("")
    rprint(
        Panel(
            f"The workflow has terminated due to an error. See the log file below for more details.",
            title = "[bold]harpy phase",
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
            title = "[bold]harpy phase",
            title_align = "left",
            border_style = "green"
            ),
        file = sys.stderr
    )

def sam_index(infile):
    """Use Samtools to index an input file, adding .bai to the end of the name"""
    if not os.path.exists(f"{infile}.bai"):
        subprocess.run(f"samtools index {infile} {infile}.bai".split())

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

rule extract_heterozygous:
    input: 
        vcf = variantfile
    output:
        outdir + "/workflow/input/vcf/{sample}.het.vcf"
    container:
        None
    message:
        "Extracting heterozygous variants: {wildcards.sample}"
    shell:
        """
        bcftools view -s {wildcards.sample} -Ou -m 2 -M 2 {input.vcf} | bcftools view -i 'GT="het"' > {output}
        """

rule extract_homozygous:
    input: 
        vcf = variantfile
    output:
        outdir + "/workflow/input/vcf/{sample}.hom.vcf"
    container:
        None
    message:
        "Extracting variants: {wildcards.sample}"
    shell:
        """
        bcftools view -s {wildcards.sample} -Ou {input.vcf} | bcftools view -i 'GT="hom"' > {output}
        """

# not the ideal way of doing this, but it works
rule index_alignments:
    input:
        bamlist
    output:
        [f"{i}.bai" for i in bamlist]
    threads:
        workflow.cores
    message:
        "Indexing alignment files"
    run:
        with multiprocessing.Pool(processes=threads) as pool:
            pool.map(sam_index, input)

if indels:
    rule genome_setup:
        input:
            genomefile
        output: 
            geno
        container:
            None
        message: 
            "Copying {input} to Genome/"
        shell: 
            """
            if (file {input} | grep -q compressed ) ;then
                # is regular gzipped, needs to be decompressed
                zcat {input} > {output}
            elif (file {input} | grep -q BGZF ); then
                # is bgzipped, also decompressed
                zcat {input} > {output}
            else
                # isn't compressed, just copied
                cp {input} {output}
            fi
            """

    rule genome_faidx:
        input: 
            geno
        output: 
            genofai
        log:
            f"Genome/{bn}.faidx.log"
        container:
            None
        message:
            "Indexing {input}"
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
        outdir + "/extractHairs/{sample}.unlinked.frags"
    log:
        outdir + "/logs/extractHairs/{sample}.unlinked.log"
    params:
        indels = indelarg,
        bx = linkarg
    conda:
        f"{envdir}/phase.yaml"
    benchmark:
        ".Benchmark/Phase/extracthairs.{sample}.txt"
    message:
        "Converting to compact fragment format: {wildcards.sample}"
    shell:
        "extractHAIRS {params} --nf 1 --bam {input.bam} --VCF {input.vcf} --out {output} > {log} 2>&1"

rule link_fragments:
    input: 
        bam       = get_alignments,
        vcf       = outdir + "/workflow/input/vcf/{sample}.het.vcf",
        fragments = outdir + "/extractHairs/{sample}.unlinked.frags"
    output:
        outdir + "/linkFragments/{sample}.linked.frags"
    log:
        outdir + "/logs/linkFragments/{sample}.linked.log"
    params:
        d = molecule_distance
    conda:
        f"{envdir}/phase.yaml"
    benchmark:
        ".Benchmark/Phase/linkfrag.{sample}.txt"
    message:
        "Linking fragments: {wildcards.sample}"
    shell:
        "LinkFragments.py --bam {input.bam} --VCF {input.vcf} --fragments {input.fragments} --out {output} -d {params} > {log} 2>&1"

rule phase_blocks:
    input:
        vcf       = outdir + "/workflow/input/vcf/{sample}.het.vcf",
        fragments = fragfile
    output: 
        blocks    = outdir + "/phaseBlocks/{sample}.blocks",
        vcf       = temp(outdir + "/phaseBlocks/{sample}.blocks.phased.VCF")
    log:
        outdir + "/logs/phaseBlocks/{sample}.blocks.phased.log"
    params: 
        prune = f"--threshold {pruning}" if pruning > 0 else "--no_prune 1",
        extra = extra
    conda:
        f"{envdir}/phase.yaml"
    benchmark:
        ".Benchmark/Phase/phase.{sample}.txt"
    message:
        "Creating phased haplotype blocks: {wildcards.sample}"
    shell:
        "HAPCUT2 --fragments {input.fragments} --vcf {input.vcf} {params} --out {output.blocks} --nf 1 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1 > {log} 2>&1"

rule compress_phaseblock:
    input:
        outdir + "/phaseBlocks/{sample}.blocks.phased.VCF"
    output:
        outdir + "/phaseBlocks/{sample}.phased.vcf.gz"
    container:
        None
    message:
        "Compressing vcf: {wildcards.sample}"
    shell:
        "bcftools view -Oz6 -o {output} --write-index {input}"

use rule compress_phaseblock as compress vcf with:
    input:
        outdir + "/workflow/input/vcf/{sample}.hom.vcf"
    output:
        outdir + "/workflow/input/gzvcf/{sample}.hom.vcf.gz"

rule merge_annotations:
    input:
        phase = outdir + "/phaseBlocks/{sample}.phased.vcf.gz",
        orig  = outdir + "/workflow/input/gzvcf/{sample}.hom.vcf.gz"
    output:
        bcf = outdir + "/annotations/{sample}.phased.annot.bcf",
        idx = outdir + "/annotations/{sample}.phased.annot.bcf.csi"
    threads:
        2
    benchmark:
        ".Benchmark/Phase/mergeAnno.{sample}.txt"
    container:
        None
    message:
        "Merging annotations: {wildcards.sample}"
    shell:
        "bcftools annotate -Ob -o {output.bcf} --write-index -a {input.phase} -c CHROM,POS,FMT/GT,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT {input.orig}"

rule merge_samples:
    input: 
        bcf = collect(outdir + "/annotations/{sample}.phased.annot.bcf", sample = samplenames),
        idx = collect(outdir + "/annotations/{sample}.phased.annot.bcf.csi", sample = samplenames)
    output:
        bcf = outdir + "/variants.phased.bcf",
        idx = outdir + "/variants.phased.bcf.csi"
    params:
        "true" if len(samplenames) > 1 else "false"
    threads:
        workflow.cores
    container:
        None
    message:
        "Combining results into {output.bcf}" if len(samplenames) > 1 else "Copying results to {output.bcf}"
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
        collect(outdir + "/phaseBlocks/{sample}.blocks", sample = samplenames)
    output:
        outdir + "/reports/blocks.summary.gz"
    params:
        outdir + "/reports/blocks.summary"
    container:
        None
    message:
        "Summarizing phasing results"
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
    conda:
        f"{envdir}/r.yaml"
    message:
        "Summarizing phasing results"
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
    message:
        "Creating record of relevant runtime parameters"
    run:
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
