bam_dir = config["seq_directory"]
genomefile = config["genomefile"]
samplenames = config["samplenames"] 
extra = config.get("extra", "") 

bn = os.path.basename(genomefile)
shell("mkdir -p Assembly")
if not os.path.exists(f"Assembly/{bn}"):
    shell(f"ln -sr {infile} Assembly/{bn}")

rule keep_validBX:
    input: bam_dir + "/{sample}.bam"
    output: "Variants/leviathan/validBX/{sample}.bx.valid.bam"
    message: "Keeping only alignments with valid BX barcodes: {wildcards.sample}"
    shell:
        """
        utilities/filterBXBAM.py --valid --input {input}
        """

rule index_alignment:
    input: "Variants/leviathan/validBX/{sample}.bx.valid.bam"
    output: "Variants/leviathan/validBX/{sample}.bx.valid.bam.bai"
    message: "Indexing barcodes: {wildcards.sample}"
    benchmark: "Benchmark/Variants/leviathan/indexbam.{sample}.txt"
    threads: 1
    shell:
        """
        sambamba index {input} {output}
        """

rule index_barcode:
    input: 
        bam = "Variants/leviathan/validBX/{sample}.bx.valid.bam",
        bai = "Variants/leviathan/validBX/{sample}.bx.valid.bam.bai"
    output: temp("Variants/leviathan/lrezIndexed/{sample}.bci")
    message: "Indexing barcodes: {wildcards.sample}"
    benchmark: "Benchmark/Variants/leviathan/indexbc.{sample}.txt"
    threads: 4
    shell:
        """
        LRez index bam -p -b {input.bam} -o {output} --threads {threads}
        """

rule leviathan_variantcall:
    input:
        bam = "Variants/leviathan/validBX/{sample}.bx.valid.bam",
        bai = "Variants/leviathan/validBX/{sample}.bx.valid.bam.bai",
        bc_idx = "Variants/leviathan/lrezIndexed/{sample}.bci",
        genome = f"Assembly/{genomefile}"
    output: vcf = pipe("Variants/leviathan/{sample}.vcf")
    log:  
        runlog = "Variants/leviathan/logs/{sample}.leviathan.log",
        candidates = "Variants/leviathan/logs/{sample}.candidates"
    message: "Calling variants: {wildcards.sample}"
    benchmark: "Benchmark/Variants/leviathan/variantcall.{sample}.txt"
    params:
        extra = extra
    threads: 3
    shell:
        """
        LEVIATHAN -b {input.bam} -i {input.bc_idx} {params} -g {input.genome} -o {output} -t {threads} --candidates {log.candidates} 2> {log.runlog}
        """

rule sort_bcf:
    input: "Variants/leviathan/{sample}.vcf"
    output: "Variants/leviathan/{sample}.bcf"
    message: "Sorting and converting to BCF: {wildcards.sample}"
    threads: 1
    params: "{wildcards.sample}"
    benchmark: "Benchmark/Variants/leviathan/sortbcf.{sample}.txt"
    shell:        
        """
        bcftools sort -Ob --output {output} {input} 2> /dev/null
        """

rule index_bcf:
    input: "Variants/leviathan/{sample}.bcf"
    output: "Variants/leviathan/{sample}.bcf.csi"
    message: "Indexing: {input}"
    benchmark: "Benchmark/Variants/leviathan/indexbcf.{sample}.txt"
    threads: 1
    shell:
        """
        bcftools index --output {output} {input}
        """

rule sv_stats:
    input: 
        bcf = "Variants/leviathan/{sample}.bcf",
        idx = "Variants/leviathan/{sample}.bcf.csi"
    output: "Variants/leviathan/stats/{sample}.sv.stats"
    message: "Getting stats for {input.bcf}"
    benchmark: "Benchmark/Variants/leviathan/stats.{sample}.txt"
    threads: 1
    shell:
        """
        bcftools stats {input.bcf} > {output}
        """

rule all_bcfs:
    input: 
        bcf = expand("Variants/leviathan/{sample}.bcf", sample = samplenames),
        index = expand("Variants/leviathan/stats/{sample}.sv.stats", sample = samplenames)
    message: "Variant calling is complete!"
    default_target: True