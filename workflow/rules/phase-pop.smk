##TODO MANUAL PRUNING OF SWITCH ERRORS
# https://github.com/vibansal/HapCUT2/blob/master/outputformat.md

bam_dir           = config["seq_directory"]
samplenames       = config["samplenames"]
variantfile       = config["variantfile"]
pruning           = config["prune"]
molecule_distance = config["molecule_distance"]
extra             = config.get("extra", "") 
outdir 			  = "Phase.noBX"if config["noBX"] else "Phase"
fragfile          = "Phase.noBX/extractHairs/{sample}.unlinked.frags" if config["noBX"] else "Phase/linkFragments/{sample}.linked.frags"
linkarg           = "--10x 0" if config["noBX"] else "--10x 1"
try:
    indelarg = "--indels 1 --ref " + config["indels"]
except:
    indelarg = ""

conda:
    os.getcwd() + "/harpyenvs/phase.yaml"

rule splitbysamplehet:
    input: 
        vcf = variantfile,
        bam = bam_dir + "/{sample}.bam"
    output:
        outdir + "/input/{sample}.het.bcf"
    message:
        "Extracting heterozygous variants: {wildcards.sample}"
    benchmark:
        ".Benchmark/Phase/splithet.{sample}.txt"
    shell:
        """
        bcftools view -s {wildcards.sample} {input.vcf} |
        awk '/^#/;/CHROM/ {{OFS="\\t"}}; !/^#/ && $10~/^0\\/1/' > {output}
        """

rule splitbysample:
    input: 
        vcf = variantfile,
        bam = bam_dir + "/{sample}.bam"
    output:
        outdir + "/input/{sample}.bcf"
    message:
        "Extracting variants: {wildcards.sample}"
    benchmark:
        ".Benchmark/Phase/split.{sample}.txt"
    shell:
        """
        bcftools view -s {wildcards.sample} {input.vcf} |
        awk '/^#/;/CHROM/ {{OFS="\\t"}}; !/^#/ &&  $10~/^0\\/0/ {{$10="0|0:"substr($10,5);print $0}}; !/^#/ && $10~/^0\\/1/; !/^#/ &&  $10~/^1\\/1/ {{$10="1|1:"substr($10,5);print $0}}; !/^#/ {{print $0}}' > {output}
        """

rule extractHairs:
    input:
        vcf = outdir + "/input/{sample}.het.bcf",
        bam = bam_dir + "/{sample}.bam"
    output:
        outdir + "/extractHairs/{sample}.unlinked.frags"
    log:
        outdir + "/extractHairs/logs/{sample}.unlinked.log"
    message:
        "Converting to compact fragment format: {wildcards.sample}"
    params:
        indels = indelarg,
        bx = linkarg
    benchmark:
        ".Benchmark/Phase/extracthairs.{sample}.txt"
    shell:
        "extractHAIRS {params} --nf 1 --bam {input.bam} --VCF {input.vcf} --out {output} 2> {log}"

rule linkFragments:
    input: 
        bam       = bam_dir + "/{sample}.bam",
        vcf       = outdir + "/input/{sample}.het.bcf",
        fragments = outdir + "/extractHairs/{sample}.unlinked.frags"
    output:
        outdir + "/linkFragments/{sample}.linked.frags"
    log:
        outdir + "/linkFragments/logs/{sample}.linked.log"
    message:
        "Linking fragments: {wildcards.sample}"
    benchmark:
        ".Benchmark/Phase/linkfrag.{sample}.txt"
    params:
        d = molecule_distance
    shell:
        "LinkFragments.py --bam {input.bam} --VCF {input.vcf} --fragments {input.fragments} --out {output} -d {params} > {log} 2>&1"

rule phaseBlocks:
    input:
        vcf       = outdir + "/input/{sample}.het.bcf",
        fragments = fragfile
    output: 
        blocks    = outdir + "/phaseBlocks/{sample}.blocks",
        vcf       = outdir + "/phaseBlocks/{sample}.blocks.phased.VCF"
    message:
        "Creating phased haplotype blocks: {wildcards.sample}"
    benchmark:
        ".Benchmark/Phase/phase.{sample}.txt"
    log:
        outdir + "/phaseBlocks/logs/{sample}.blocks.phased.log"
    params: 
        prune = f"--threshold {pruning}" if pruning > 0 else "--no_prune 1",
        extra = extra
    shell:
        "HAPCUT2 --fragments {input.fragments} --vcf {input.vcf} {params} --out {output.blocks} --nf 1 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1 > {log} 2>&1"

rule createAnnotations:
    input:
        outdir + "/phaseBlocks/{sample}.blocks.phased.VCF"
    output:
        outdir + "/annotations/{sample}.annot.gz"
    message:
        "Creating annotation files: {wildcards.sample}"
    benchmark:
        ".Benchmark/Phase/createAnno.{sample}.txt"
    shell:
        "bcftools query -f \"%CHROM\\t%POS[\\t%GT\\t%PS\\t%PQ\\t%PD]\\n\" {input} | bgzip -c > {output}"

rule indexAnnotations:
    input:
        outdir + "/annotations/{sample}.annot.gz"
    output:
        outdir + "/annotations/{sample}.annot.gz.tbi"
    message:
        "Indexing {wildcards.sample}.annot.gz"
    benchmark:
        ".Benchmark/Phase/indexAnno.{sample}.txt"
    shell: 
        "tabix -b 2 -e 2 {input}"

rule headerfile:
    output:
        outdir + "/input/header.names"
    message:
        "Creating additional header file"
    benchmark:
        ".Benchmark/Phase/headerfile.txt"
    run:
        with open(output[0], "w") as fout:
            _ = fout.write('##INFO=<ID=HAPCUT,Number=0,Type=Flag,Description="The haplotype was created with Hapcut2">\n')
            _ = fout.write('##FORMAT=<ID=GX,Number=1,Type=String,Description="Haplotype">\n')
            _ = fout.write('##FORMAT=<ID=PS,Number=1,Type=Integer,Description="ID of Phase Set for Variant">\n')
            _ = fout.write('##FORMAT=<ID=PQ,Number=1,Type=Integer,Description="Phred QV indicating probability that this variant is incorrectly phased relative to the haplotype">\n')
            _ = fout.write('##FORMAT=<ID=PD,Number=1,Type=Integer,Description="phased Read Depth">')

rule mergeAnnotations:
    input:
        annot   = outdir + "/annotations/{sample}.annot.gz",
        idx     = outdir + "/annotations/{sample}.annot.gz.tbi",
        orig    = outdir + "/input/{sample}.bcf",
        headers = outdir + "/input/header.names"
    output:
        bcf = outdir + "/annotations_merge/{sample}.phased.annot.bcf",
        idx = outdir + "/annotations_merge/{sample}.phased.annot.bcf.csi"

    message:
        "Merging annotations: {wildcards.sample}"
    threads:
        2
    benchmark:
        ".Benchmark/Phase/mergeAnno.{sample}.txt"
    shell:
        """
        bcftools annotate -h {input.headers} -a {input.annot} {input.orig} -c CHROM,POS,FMT/GX,FMT/PS,FMT/PQ,FMT/PD -m +HAPCUT |
        awk '!/<ID=GX/' |
        sed 's/:GX:/:GT:/' |
        #bcftools view -Ob --write-index -o {output.bcf} -
        bcftools view -Ob -o {output.bcf} -
        bcftools index {output.bcf}
        """

rule mergeSamples:
    input: 
        bcf = expand(outdir + "/annotations_merge/{sample}.phased.annot.bcf", sample = samplenames),
        idx = expand(outdir + "/annotations_merge/{sample}.phased.annot.bcf.csi", sample = samplenames)
    output:
        bcf = outdir + "/variants.phased.bcf",
        idx = outdir + "/variants.phased.bcf.csi"
    message:
        "Combining samples into a single BCF file"
    benchmark:
        ".Benchmark/Phase/mergesamples.txt"
    threads:
        30
    shell:
        """
        bcftools merge --threads {threads} -Ob {input.bcf} > {output.bcf}
        bcftools index {output.bcf}
        """
        #"bcftools merge --threads {threads} --Ob --write-index {input.bcf} > {output}"

rule summarize_blocks:
    input:
        expand(outdir + "/phaseBlocks/{sample}.blocks", sample = samplenames)
    output:
        outdir + "/reports/blocks.summary.gz"
    message:
        "Summarizing phasing results"
    params:
        outdir + "/reports/blocks.summary"
    shell:
        """
        echo -e "sample\\tcontig\\tn_snp\\tpos_start\\tblock_length" > {params}
        for i in {input}; do
            parsePhaseBlocks.py -i $i >> {params}
        done
        gzip {params}
        """

rule phase_report:
    input:
        outdir + "/reports/blocks.summary.gz"
    output:
        outdir + "/reports/phase.html"
    conda:
        os.getcwd() + "/harpyenvs/r-env.yaml"
    message:
        "Summarizing phasing results"
    script:
        "reportHapCut2.Rmd"

rule log_runtime:
    output:
        outdir + "/logs/harpy.phase.log"
    message:
        "Creating record of relevant runtime parameters"
    params:
        prune = f"--threshold {pruning}" if pruning > 0 else "--no_prune 1",
        extra = extra
    run:
        with open(output[0], "w") as f:
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
            _ = f.write("    awk '!/<ID=GX/' |\n")
            _ = f.write("    sed 's/:GX:/:GT:/' |\n")
            _ = f.write("    bcftools view -Ob -o sample.annot.bcf -\n")
            _ = f.write("    bcftools merge --output-type b samples.annot.bcf\n\n")
            _ = f.write("The header.file of extra vcf tags:\n")
            _ = f.write('    ##INFO=<ID=HAPCUT,Number=0,Type=Flag,Description="The haplotype was created with Hapcut2">\n')
            _ = f.write('    ##FORMAT=<ID=GX,Number=1,Type=String,Description="Haplotype">\n')
            _ = f.write('    ##FORMAT=<ID=PS,Number=1,Type=Integer,Description="ID of Phase Set for Variant">\n')
            _ = f.write('    ##FORMAT=<ID=PQ,Number=1,Type=Integer,Description="Phred QV indicating probability that this variant is incorrectly phased relative to the haplotype">\n')
            _ = f.write('    ##FORMAT=<ID=PD,Number=1,Type=Integer,Description="phased Read Depth">\n')


rule indexFinal:
    default_target: True
    input:
        outdir + "/variants.phased.bcf",
        outdir + "/logs/harpy.phase.log",
        outdir + "/reports/phase.html"
    benchmark:
        ".Benchmark/Phase/finalindex.txt"
    message:
        "Phasing is complete!"

