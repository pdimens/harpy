bam_dir           = config["seq_directory"]
samplenames       = config["samplenames"]
variantfile       = config["variantfile"]
pruning           = config["prune"]
molecule_distance = config["molecule_distance"]
extra             = config.get("extra", "") 

rule splitbysamplehet:
    input: 
        vcf = variantfile,
        bam = bam_dir + "/{sample}.bam"
    output:
        "Phase/input/{sample}.het.bcf"
    message:
        "Extracting heterozygous variants: {wildcards.sample}"
    benchmark:
        "Benchmark/Phase/splithet.{sample}.txt"
    threads: 1
    shell:
        """
        {
            bcftools view -s {wildcards.sample} -i 'INFO/INFO_SCORE >= 0.2' {input.vcf} 
        } || {
            bcftools view -s {wildcards.sample} {input.vcf} 
        } |
        awk '/^#/;/CHROM/ {{OFS="\\t"}}; !/^#/ && $10~/^0\\/1/' > {output}
        """

rule splitbysample:
    input: 
        vcf = variantfile,
        bam = bam_dir + "/{sample}.bam"
    output:
        "Phase/input/{sample}.bcf"
    message:
        "Extracting variants: {wildcards.sample}"
    benchmark:
        "Benchmark/Phase/split.{sample}.txt"
    threads: 1
    shell:
        """
        {
            bcftools view -s {wildcards.sample} -i 'INFO/INFO_SCORE >= 0.2' {input.vcf} |
        } || {
            bcftools view -s {wildcards.sample} {input.vcf} |
        } | 
        awk '/^#/;/CHROM/ {{OFS="\\t"}}; !/^#/ &&  $10~/^0\\/0/ {{$10="0|0:"substr($10,5);print $0}}; !/^#/ && $10~/^0\\/1/; !/^#/ &&  $10~/^1\\/1/ {{$10="1|1:"substr($10,5);print $0}}; !/^#/ {{print $0}}' > {output}
        """

rule extractHairs:
    input:
        vcf = "Phase/input/{sample}.het.bcf",
        bam = bam_dir + "/{sample}.bam"
    output:
        "Phase/extractHairs/{sample}.unlinked.frags"
    log:
        "Phase/extractHairs/logs/{sample}.unlinked.log"
    message:
        "Converting to compact fragment format: {wildcards.sample}"
    benchmark:
        "Benchmark/Phase/extracthairs.{sample}.txt"
    threads: 1
    shell:
        "extractHAIRS --10X 1 --nf 1 --bam {input.bam} --VCF {input.vcf} --out {output} 2> {log}"

rule linkFragments:
    input: 
        bam       = bam_dir + "/{sample}.bam",
        vcf       = "Phase/input/{sample}.het.bcf",
        fragments = "Phase/extractHairs/{sample}.unlinked.frags"
    output:
        "Phase/linkFragments/{sample}.linked.frags"
    log:
        "Phase/linkFragments/logs/{sample}.linked.log"
    message:
        "Linking fragments: {wildcards.sample}"
    benchmark:
        "Benchmark/Phase/linkfrag.{sample}.txt"
    params:
        d = molecule_distance
    shell:
        "LinkFragments.py  --bam {input.bam} --VCF {input.vcf} --fragments {input.fragments} --out {output} -d {params} > {log} 2>&1"

rule phaseBlocks:
    input:
        vcf       = "Phase/input/{sample}.het.bcf",
        fragments = "Phase/linkFragments/{sample}.linked.frags"
    output: 
        blocks    = "Phase/phaseBlocks/{sample}.blocks",
        vcf       = "Phase/phaseBlocks/{sample}.blocks.phased.VCF"
    message:
        "Creating phased haplotype blocks: {wildcards.sample}"
    benchmark:
        "Benchmark/Phase/phase.{sample}.txt"
    log:
        "Phase/phaseBlocks/logs/{sample}.blocks.phased.log"
    params: 
        prune = f"--threshold {pruning}" if pruning > 0 else "--no_prune 1",
        extra = extra
    threads: 1
    shell:
        "HAPCUT2 --fragments {input.fragments} --vcf {input.vcf} {params} --out {output.blocks} --nf 1 {params} --error_analysis_mode 1 --call_homozygous 1 --outvcf 1 2> {log}"

rule createAnnotations:
    input:
        "Phase/phaseBlocks/{sample}.blocks.phased.VCF"
    output:
        "Phase/annotations/{sample}.annot.gz"
    message:
        "Creating annotation files: {wildcards.sample}"
    benchmark:
        "Benchmark/Phase/createAnno.{sample}.txt"
    shell:
        "bcftools query -f \"%CHROM\\t%POS[\\t%GT\\t%PS\\t%PQ\\t%PD]\\n\" {input} | bgzip -c > {output}"

rule indexAnnotations:
    input:
        "Phase/annotations/{sample}.annot.gz"
    output:
        "Phase/annotations/{sample}.annot.gz.tbi"
    message:
        "Indexing {wildcards.sample}.annot.gz"
    benchmark:
        "Benchmark/Phase/indexAnno.{sample}.txt"
    shell: 
        "tabix -b 2 -e 2 {input}"

rule headerfile:
    output:
        "Phase/input/header.names"
    message:
        "Creating additional header file"
    benchmark:
        "Benchmark/Phase/headerfile.txt"
    run:
        with open(output[0], "w") as fout:
            fout.write('##INFO=<ID=HAPCUT,Number=0,Type=Flag,Description="The haplotype was created with Hapcut2">\n')
            fout.write('##FORMAT=<ID=GX,Number=1,Type=String,Description="Haplotype">\n')
            fout.write('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="ID of Phase Set for Variant">\n')
            fout.write('##FORMAT=<ID=PQ,Number=1,Type=Integer,Description="Phred QV indicating probability that this variant is incorrectly phased relative to the haplotype">\n')
            fout.write('##FORMAT=<ID=PD,Number=1,Type=Integer,Description="phased Read Depth">')

rule mergeAnnotations:
    input:
        annot   = "Phase/annotations/{sample}.annot.gz",
        idx     = "Phase/annotations/{sample}.annot.gz.tbi",
        orig    = "Phase/input/{sample}.bcf",
        headers = "Phase/input/header.names"
    output:
        "Phase/annotations_merge/{sample}.phased.annot.bcf"
    message:
        "Merging annotations: {wildcards.sample}"
    benchmark:
        "Benchmark/Phase/mergeAnno.{sample}.txt"
    shell:
        "bcftools annotate -h {input.headers} -a {input.annot} {input.orig} -c CHROM,POS,FMT/GX,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT |  awk '!/<ID=GX/' | sed 's/:GX:/:GT:/' | bcftools view -Ob -o {output} -"
        
rule indexAnnotations2:
    input:
        "Phase/annotations_merge/{sample}.phased.annot.bcf"
    output:
        "Phase/annotations_merge/{sample}.phased.annot.bcf.csi"
    message:
        "Indexing annotations: {wildcards.sample}"
    benchmark:
        "Benchmark/Phase/indexAnno.{sample}.txt"
    shell:
        "bcftools index {input}"

rule mergeSamples:
    input: 
        bcf = expand("Phase/annotations_merge/{sample}.phased.annot.bcf", sample = samplenames),
        idx = expand("Phase/annotations_merge/{sample}.phased.annot.bcf.csi", sample = samplenames)
    output:
        "Phase/variants.phased.bcf"
    message:
        "Combinging samples into a single BCF file"
    benchmark:
        "Benchmark/Phase/mergesamples.txt"
    threads: 30
    shell:
        "bcftools merge --threads {threads} --output-type b {input.bcf} > {output}"

rule indexFinal:
    input:
        "Phase/variants.phased.bcf"
    output:
        "Phase/variants.phased.bcf.csi"
    benchmark:
        "Benchmark/Phase/finalindex.txt"
    message:
        "Phasing is complete!"
    default_target: True
    shell: 
        "bcftools index {input}"
