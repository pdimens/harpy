from rich import print as rprint
from rich.panel import Panel
import sys

bam_dir           = config["seq_directory"]
samplenames       = config["samplenames"]
variantfile       = config["variantfile"]
pruning           = config["prune"]
molecule_distance = config["molecule_distance"]
extra             = config.get("extra", "") 
linkarg           = "--10x 0" if config["noBX"] else "--10x 1"
outdir 			  = "Phase.noBX"if config["noBX"] else "Phase"
fragfile          = "Phase.noBX/extractHairs/{sample}.unlinked.frags" if config["noBX"] else "Phase/linkFragments/{sample}.linked.frags"
skipreports = config["skipreports"]

try:
    indelarg = "--indels 1 --ref " + config["indels"]
except:
    indelarg = ""

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


rule splitbysample:
    input: 
        vcf = variantfile,
        bam = bam_dir + "/{sample}.bam"
    output:
        temp(outdir + "/input/{sample}.bcf")
    benchmark:
        ".Benchmark/Phase/split.{sample}.txt"
    wildcard_constraints:
        sample = "[a-zA-Z0-9_-]*"
    message:
        "Extracting variants: {wildcards.sample}"
    shell:
        """
        bcftools view -s {wildcards.sample} {input.vcf} |
        awk '/^#/;/CHROM/ {{OFS="\\t"}}; !/^#/ && $10~/^0\\/1/' > {output}
        """

rule extractHairs:
    input:
        vcf = "Phase/input/{sample}.bcf",
        bam = bam_dir + "/{sample}.bam"
    output:
        outdir + "/extractHairs/{sample}.unlinked.frags"
    log:
        outdir + "/extractHairs/logs/{sample}.unlinked.log"
    params:
        linkarg
    wildcard_constraints:
        sample = "[a-zA-Z0-9_-]*"
    benchmark:
        ".Benchmark/Phase/extracthairs.{sample}.txt"
    conda:
        os.getcwd() + "/.harpy_envs/phase.yaml"
    message:
        "Converting to compact fragment format: {wildcards.sample}"
    shell:
        "extractHAIRS {params} --nf 1 --bam {input.bam} --VCF {input.vcf} --out {output} 2> {log}"

rule linkFragments:
    input: 
        bam       = bam_dir + "/{sample}.bam",
        vcf       = outdir + "/input/{sample}.bcf",
        fragments = outdir + "/extractHairs/{sample}.unlinked.frags"
    output:
        outdir + "/linkFragments/{sample}.linked.frags"
    log:
        outdir + "/linkFragments/logs/{sample}.linked.log"
    benchmark:
        ".Benchmark/Phase/linkfrag.{sample}.txt"
    params:
        d = molecule_distance
    conda:
        os.getcwd() + "/.harpy_envs/phase.yaml"
    message:
        "Linking fragments: {wildcards.sample}"
    shell:
        "LinkFragments.py --bam {input.bam} --VCF {input.vcf} --fragments {input.fragments} --out {output} -d {params} > {log} 2>&1"


rule phaseBlocks:
    input:
        vcf       = outdir + "/input/{sample}.bcf",
        fragments = fragfile
    output: 
        blocks    = outdir + "/phaseBlocks/{sample}.blocks",
        vcf       = temp(outdir + "/phaseBlocks/{sample}.blocks.phased.VCF")
    benchmark:
        ".Benchmark/Phase/phase.{sample}.txt"
    log:
        outdir + "/phaseBlocks/logs/{sample}.blocks.phased.log"
    params: 
        prune = f"--threshold {pruning}" if pruning > 0 else "--no_prune 1",
        extra = extra
    conda:
        os.getcwd() + "/.harpy_envs/phase.yaml"
    message:
        "Creating phased haplotype blocks: {wildcards.sample}"
    shell:
        "HAPCUT2 --fragments {input.fragments} --vcf {input.vcf} {params} --out {output.blocks} --nf 1 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1 2> {log}"

rule compress_vcf:
    input:
        vcf = outdir + "/phaseBlocks/{sample}.blocks.phased.VCF"
    output: 
        idx = outdir + "/phaseBlocks/{sample}.blocks.phased.VCF.gz",
    message:
        "Compressing: {wildcards.sample}"
    shell:
        "bgzip {input}"

rule index_vcf:
    input:
        vcf = outdir + "/phaseBlocks/{sample}.blocks.phased.VCF.gz"
    output: 
        idx = outdir + "/phaseBlocks/{sample}.blocks.phased.VCF.gz.tbi",
    message:
        "Compressing: {wildcards.sample}"
    shell:
        "tabix -b 2 -e 2 {input}"

rule mergeSamples:
    input: 
        vcf = expand(outdir + "/phaseBlocks/{sample}.blocks.phased.VCF.gz", sample = samplenames),
        idx = expand(outdir + "/phaseBlocks/{sample}.blocks.phased.VCF.gz.tbi", sample = samplenames)
    output:
        outdir + "/variants.phased.bcf"
    benchmark:
        ".Benchmark/Phase/mergesamples.txt"
    threads:
        30
    message:
        "Combining samples into a single BCF file"
    shell:
        """
        bcftools merge --threads {threads} --output-type b --write-index -o {output} {input.vcf}
        """

rule merge_blocks:
    input:
        expand(outdir + "/phaseBlocks/{sample}.blocks", sample = samplenames)
    output:
        outdir + "/phased.blocks.gz"
    message:
        "Summarizing phasing results"
    shell:
        "cat {input} | gzip > {output}"

rule phase_report:
    input:
        outdir + "/phased.blocks.gz"
    output:
        outdir + "/reports/phase.html"
    conda:
        os.getcwd() + "/.harpy_envs/r-env.yaml"
    message:
        "Summarizing phasing results"
    script:
        "report/HapCut2.Rmd"


rule log_workflow:
    default_target: True
    input:
        vcf = outdir + "/variants.phased.bcf",
        report = outdir + "/reports/phase.html" if not skipreports else []
    output:
        outdir + "/workflow/phase.workflow.summary"
    params:
        prune = f"--threshold {pruning} " if pruning > 0 else "--no_prune 1 ",
        extra = extra
    message:
        "Summarizing the workflow: {output}"
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy phase module ran using these parameters:\n\n")
            _ = f.write(f"The provided variant file: {variantfile}\n")
            _ = f.write(f"The directory with alignments: {bam_dir}\n")
            _ = f.write("The variant file was split by sample and preprocessed using:\n")
            _ = f.write("""    bcftools view -s SAMPLE | awk '/^#/;/CHROM/ OFS="\\t"; !/^#/ && $10~/^0\\/1/'\n\n""")
            _ = f.write("Phasing was performed using the components of HapCut2:\n")
            _ = f.write(f"    extractHAIRS {linkarg} --nf 1 --bam sample.bam --VCF sample.vcf --out sample.unlinked.frags\n")
            _ = f.write(f"    LinkFragments.py --bam sample.BAM --VCF sample.vcf --fragments sample.unlinked.frags --out sample.linked.frags -d {molecule_distance}\n")
            _ = f.write(f"    HAPCUT2 --fragments sample.linked.frags --vcf sample.vcf --out sample.blocks --nf 1 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1 {params[0]} {params[1]}\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")