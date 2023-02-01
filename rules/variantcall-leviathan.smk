bam_dir = config["seq_directory"]
genomefile = config["genomefile"]
samplenames = config["samplenames"] 
extra = config["extra"]

rule index_alignment:
    input: bam_dir + "/{sample}.bam"
    output: bam_dir + "/{sample}.bam.bai"
    message: "Indexing barcodes: {wildcards.sample}"
    threads: 1
    shell:
        """
        sambamba index {input} {output}
        """

rule index_barcode:
    input: 
        bam = bam_dir + "/{sample}.bam",
        bai = bam_dir + "/{sample}.bam.bai"
    output: temp("Variants/leviathan/lrezIndexed/{sample}.bci")
    message: "Indexing barcodes: {wildcards.sample}"
    threads: 4
    shell:
        """
        LRez index bam -p -b {input} -o {output}
        """

rule leviathan_variantcall:
    input:
        bam = bam_dir + "/{sample}" + ".bam",
        bai = bam_dir + "/{sample}" + ".bam.bai",
        bc_idx = "Variants/leviathan/lrezIndexed/{sample}.bci",
        genome = genomefile
    output: vcf = pipe("Variants/leviathan/{sample}.vcf")
    log:  
        runlog = "Variants/leviathan/logs/{sample}.leviathan.log",
        candidates = "Variants/leviathan/logs/{sample}.candidates"
    message: "Calling variants: {wildcards.sample}"
    params:
        extra = extra
    threads: 3
    shell:
        """
        LEVIATHAN -b {input.bam} -i {input.bc_idx} {params} -g {input.genome} -o {output} -t {threads} --candidates {log.candidates} 2> {log.runlog}
        """

rule convert_bcf:
    input: "Variants/leviathan/{sample}.vcf"
    output: temp("Variants/leviathan/{sample}.bcf")
    message: "Covnerting to BCF: {input}"
    threads: 1
    params: "{wildcards.sample}"
    shell:        
        """
        bcftools sort -Ob --output {output} {input} 2> /dev/null
        """

rule index_bcf:
    input: "Variants/leviathan/{sample}.bcf"
    output: "Variants/leviathan/{sample}.bcf.csi"
    message: "Indexing: {input}"
    threads: 1
    shell:
        """
        bcftools index --output {output} {input}
        """

rule sv_stats:
    input: 
        bcf = "Variants/leviathan/{sample}.bcf",
        idx = "Variants/leviathan/{sample}.bcf.csi",
        genome = genomefile
    output: "Variants/leviathan/stats/{sample}.sv.stats"
    message: "Getting stats for {input.bcf}"
    threads: 1
    shell:
        """
        bcftools stats --fasta-ref {input.genome} {input.bcf} > {output}
        """

rule all_bcfs:
    input: 
        bcf = expand("Variants/leviathan/{sample}.bcf", sample = samplenames),
        index = expand("Variants/leviathan/stats/{sample}.sv.stats", sample = samplenames)
    message: "Variant calling is complete!"
    default_target: True