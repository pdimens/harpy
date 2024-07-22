containerized: "docker://pdimens/harpy:latest"

from rich import print as rprint
from rich.panel import Panel
import sys

##TODO MANUAL PRUNING OF SWITCH ERRORS
# https://github.com/vibansal/HapCUT2/blob/master/outputformat.md

bam_dir           = config["seq_directory"]
samplenames       = config["samplenames"]
variantfile       = config["variantfile"]
pruning           = config["prune"]
molecule_distance = config["molecule_distance"]
extra             = config.get("extra", "") 
outdir 			  = config["output_directory"]
envdir      = os.getcwd() + "/.harpy_envs"

if config["noBX"]:
    fragfile = outdir + "/extractHairs/{sample}.unlinked.frags"
else:
    fragfile =  outdir + "/linkFragments/{sample}.linked.frags"

linkarg     = "--10x 0" if config["noBX"] else "--10x 1"
skipreports = config["skip_reports"]

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

rule split_by_samplehet:
    input: 
        vcf = variantfile,
        bam = bam_dir + "/{sample}.bam"
    output:
        outdir + "/workflow/input/vcf/{sample}.het.vcf"
    benchmark:
        ".Benchmark/Phase/splithet.{sample}.txt"
    container:
        None
    message:
        "Extracting heterozygous variants: {wildcards.sample}"
    shell:
        """
        bcftools view -s {wildcards.sample} {input.vcf} |
        awk '/^#/;/CHROM/ {{OFS="\\t"}}; !/^#/ && $10~/^0\\/1/' > {output}
        """

rule split_by_sample:
    input: 
        vcf = variantfile,
        bam = bam_dir + "/{sample}.bam"
    output:
        outdir + "/workflow/input/vcf/{sample}.vcf"
    benchmark:
        ".Benchmark/Phase/split.{sample}.txt"
    container:
        None
    message:
        "Extracting variants: {wildcards.sample}"
    shell:
        """
        bcftools view -s {wildcards.sample} {input.vcf} |
        awk '/^#/;/CHROM/ {{OFS="\\t"}}; !/^#/ &&  $10~/^0\\/0/ {{$10="0|0:"substr($10,5);print $0}}; !/^#/ && $10~/^0\\/1/; !/^#/ &&  $10~/^1\\/1/ {{$10="1|1:"substr($10,5);print $0}}; !/^#/ {{print $0}}' > {output}
        """

rule index_alignment:
    input:
        bam_dir + "/{sample}.bam"
    output:
        bam_dir + "/{sample}.bam.bai"
    container:
        None
    message:
        "Indexing alignment: {wildcards.sample}"
    shell:
        "samtools index {input} {output} 2> /dev/null"

rule extract_hairs:
    input:
        vcf = outdir + "/workflow/input/vcf/{sample}.het.vcf",
        bam = bam_dir + "/{sample}.bam",
        bai = bam_dir + "/{sample}.bam.bai"
    output:
        outdir + "/extractHairs/{sample}.unlinked.frags"
    log:
        outdir + "/extractHairs/logs/{sample}.unlinked.log"
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
        "extractHAIRS {params} --nf 1 --bam {input.bam} --VCF {input.vcf} --out {output} 2> {log}"

rule link_fragments:
    input: 
        bam       = bam_dir + "/{sample}.bam",
        vcf       = outdir + "/workflow/input/vcf/{sample}.het.vcf",
        fragments = outdir + "/extractHairs/{sample}.unlinked.frags"
    output:
        outdir + "/linkFragments/{sample}.linked.frags"
    log:
        outdir + "/linkFragments/logs/{sample}.linked.log"
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
        vcf       = outdir + "/phaseBlocks/{sample}.blocks.phased.VCF"
    log:
        outdir + "/phaseBlocks/logs/{sample}.blocks.phased.log"
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

rule create_annotations:
    input:
        outdir + "/phaseBlocks/{sample}.blocks.phased.VCF"
    output:
        outdir + "/annotations/{sample}.annot.gz"
    container:
        None
    message:
        "Creating annotation files: {wildcards.sample}"
    shell:
        "bcftools query -f \"%CHROM\\t%POS[\\t%GT\\t%PS\\t%PQ\\t%PD]\\n\" {input} | bgzip -c > {output}"

rule index_annotations:
    input:
        outdir + "/annotations/{sample}.annot.gz"
    output:
        outdir + "/annotations/{sample}.annot.gz.tbi"
    container:
        None
    message:
        "Indexing {wildcards.sample}.annot.gz"
    shell: 
        "tabix -b 2 -e 2 {input}"

rule headerfile:
    output:
        outdir + "/workflow/input/header.names"
    message:
        "Creating additional header file"
    run:
        with open(output[0], "w") as fout:
            _ = fout.write('##INFO=<ID=HAPCUT,Number=0,Type=Flag,Description="The haplotype was created with Hapcut2">\n')
            _ = fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype/Haplotype">\n')
            _ = fout.write('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="ID of Phase Set for Variant">\n')
            _ = fout.write('##FORMAT=<ID=PQ,Number=1,Type=Integer,Description="Phred QV indicating probability that this variant is incorrectly phased relative to the haplotype">\n')
            _ = fout.write('##FORMAT=<ID=PD,Number=1,Type=Integer,Description="phased Read Depth">')

rule merge_annotations:
    input:
        annot   = outdir + "/annotations/{sample}.annot.gz",
        idx     = outdir + "/annotations/{sample}.annot.gz.tbi",
        orig    = outdir + "/workflow/input/vcf/{sample}.vcf",
        headers = outdir + "/workflow/input/header.names"
    output:
        bcf = outdir + "/annotations_merge/{sample}.phased.annot.bcf",
        idx = outdir + "/annotations_merge/{sample}.phased.annot.bcf.csi"
    threads:
        2
    benchmark:
        ".Benchmark/Phase/mergeAnno.{sample}.txt"
    container:
        None
    message:
        "Merging annotations: {wildcards.sample}"
    shell:
        """
        bcftools annotate -h {input.headers} -a {input.annot} {input.orig} -c CHROM,POS,FMT/GT,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT |
            awk '!/<ID=GX/' |
            sed 's/:GX:/:GT:/' |
            bcftools view -Ob --write-index -o {output.bcf} -
        """

rule merge_samples:
    input: 
        bcf = collect(outdir + "/annotations_merge/{sample}.phased.annot.bcf", sample = samplenames),
        idx = collect(outdir + "/annotations_merge/{sample}.phased.annot.bcf.csi", sample = samplenames)
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
            _ = f.write("The harpy phase module ran using these parameters:\n\n")
            _ = f.write(f"The provided variant file: {variantfile}\n")
            _ = f.write(f"The directory with alignments: {bam_dir}\n")
            _ = f.write("The variant file was split by sample and preprocessed using:\n")
            _ = f.write("""    bcftools view -s SAMPLE | awk '/^#/;/CHROM/ OFS="\\t"; !/^#/ && $10~/^0\\/1/'\n\n""")
            _ = f.write("Phasing was performed using the components of HapCut2:\n")
            _ = f.write("    extractHAIRS " + linkarg + " --nf 1 --bam sample.bam --VCF sample.vcf --out sample.unlinked.frags\n")
            _ = f.write("    LinkFragments.py --bam sample.bam --VCF sample.vcf --fragments sample.unlinked.frags --out sample.linked.frags -d " + f"{molecule_distance}" + "\n")
            _ = f.write("    HAPCUT2 --fragments sample.linked.frags --vcf sample.vcf --out sample.blocks --nf 1 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1" + f" {params[0]} {params[1]}" + "\n\n")
            _ = f.write("Variant annotation was performed using:\n")
            _ = f.write("    bcftools query -f \"%CHROM\\t%POS[\\t%GT\\t%PS\\t%PQ\\t%PD]\\n\" sample.vcf | bgzip -c\n")
            _ = f.write("    bcftools annotate -h header.file -a sample.annot sample.bcf -c CHROM,POS,FMT/GX,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT |\n")
            _ = f.write("        awk '!/<ID=GX/' |\n")
            _ = f.write("        sed 's/:GX:/:GT:/' |\n")
            _ = f.write("        bcftools view -Ob -o sample.annot.bcf -\n")
            _ = f.write("    bcftools merge --output-type b samples.annot.bcf\n\n")
            _ = f.write("The header.file of extra vcf tags:\n")
            _ = f.write('    ##INFO=<ID=HAPCUT,Number=0,Type=Flag,Description="The haplotype was created with Hapcut2">\n')
            _ = f.write('    ##FORMAT=<ID=GX,Number=1,Type=String,Description="Haplotype">\n')
            _ = f.write('    ##FORMAT=<ID=PS,Number=1,Type=Integer,Description="ID of Phase Set for Variant">\n')
            _ = f.write('    ##FORMAT=<ID=PQ,Number=1,Type=Integer,Description="Phred QV indicating probability that this variant is incorrectly phased relative to the haplotype">\n')
            _ = f.write('    ##FORMAT=<ID=PD,Number=1,Type=Integer,Description="phased Read Depth">\n')
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")