bam_dir           = config["seq_directory"]
samplenames       = config["samplenames"]
variantfile       = config["variantfile"]
pruning           = config["prune"]
molecule_distance = config["molecule_distance"]
extra             = config.get("extra", "") 
linkarg           = "--10x 0" if config["noBX"] else "--10x 1"
outdir 			  = "Phase.noBX"if config["noBX"] else "Phase"
fragfile          = "Phase.noBX/extractHairs/{sample}.unlinked.frags" if config["noBX"] else "Phase/linkFragments/{sample}.linked.frags"

rule splitbysample:
    input: 
        vcf = variantfile,
        bam = bam_dir + "/{sample}.bam"
    output:
        temp(outdir + "/input/{sample}.bcf")
    message:
        "Extracting variants: {wildcards.sample}"
    benchmark:
        "Benchmark/Phase/split.{sample}.txt"
    wildcard_constraints:
        sample = "[a-zA-Z0-9_-]*"
    threads: 1
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
    message:
        "Converting to compact fragment format: {wildcards.sample}"
    params:
        linkarg
    wildcard_constraints:
        sample = "[a-zA-Z0-9_-]*"
    benchmark:
        "Benchmark/Phase/extracthairs.{sample}.txt"
    threads: 1
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
    message:
        "Linking fragments: {wildcards.sample}"
    benchmark:
        "Benchmark/Phase/linkfrag.{sample}.txt"
    params:
        d = molecule_distance
    shell:
        "LinkFragments.py --bam {input.bam} --VCF {input.vcf} --fragments {input.fragments} --out {output} -d {params} > {log} 2>&1"


rule phaseBlocks:
    input:
        vcf       = outdir + "/input/{sample}.bcf",
        fragments = fragfile
    output: 
        blocks    = outdir + "/phaseBlocks/{sample}.blocks",
        vcf       = temp(outdir + "/phaseBlocks/{sample}.blocks.phased.VCF")
    message:
        "Creating phased haplotype blocks: {wildcards.sample}"
    benchmark:
        "Benchmark/Phase/phase.{sample}.txt"
    log:
        outdir + "/phaseBlocks/logs/{sample}.blocks.phased.log"
    params: 
        prune = f"--threshold {pruning}" if pruning > 0 else "--no_prune 1",
        extra = extra
    threads: 1
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
    message:
        "Combinging samples into a single BCF file"
    benchmark:
        "Benchmark/Phase/mergesamples.txt"
    threads:
        30
    shell:
        """
        #bcftools merge --threads {threads} --output-type b --write-index {input.vcf} > {output}
        bcftools merge --threads {threads} --output-type b {input.vcf} > {output}
        """

rule merge_blocks:
    input:
        expand(outdir + "/phaseBlocks/{sample}.blocks", sample = samplenames)
    output:
        outdir + "/phased.blocks.gz"
    message:
        "Summarizing phasing results"
    shell:
        "cat {input} | gzip > output"

rule phase_report:
    input:
        outdir + "/phased.blocks.gz"
    output:
        outdir + "/reports/phase.html"
    message:
        "Summarizing phasing results"
    script:
        "reportHapCut2.Rmd"

rule log_runtime:
    output:
        outdir + "/logs/harpy.phase.log"
    message:
        "Creating record of relevant runtime parameters: {output}"
    params:
        links = linkarg,
        d =  molecule_distance,
        prune = f"--threshold {pruning} " if pruning > 0 else "--no_prune 1 ",
        extra = extra
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy phase module ran using these parameters:\n\n")
            _ = f.write(f"The provided variant file: {variantfile}\n")
            _ = f.write(f"The directory with alignments: {bam_dir}\n")
            _ = f.write("The variant file was split by sample and preprocessed using:\n")
            _ = f.write("""\tbcftools view -s SAMPLE | awk '/^#/;/CHROM/ OFS="\\t"; !/^#/ && $10~/^0\\/1/'\n\n""")
            _ = f.write("Phasing was performed using the components of HapCut2:\n")
            _ = f.write(f"\textractHAIRS {params[0]} --nf 1 --bam sample.bam --VCF sample.vcf --out sample.unlinked.frags\n")
            _ = f.write(f"\tLinkFragments.py --bam sample.BAM --VCF sample.vcf --fragments sample.unlinked.frags --out sample.linked.frags -d {params[1]}\n")
            _ = f.write(f"\tHAPCUT2 --fragments sample.linked.frags --vcf sample.vcf --out sample.blocks --nf 1 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1 {params[2]} {params[3]}\n")

rule indexFinal:
    input:
        bcf      = outdir + "/variants.phased.bcf",
        runlog   = outdir + "/logs/harpy.phase.log",
        phaserpt = outdir + "/reports/phase.html"
    output:
        outdir + "/variants.phased.bcf.csi"
    benchmark:
        "Benchmark/Phase/finalindex.txt"
    message:
        "Phasing is complete!"
    default_target: True
    shell: 
        "bcftools index {input.bcf}"