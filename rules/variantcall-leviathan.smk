# user specified configs
bam_dir = config["seq_directory"]
genomefile = config["genome_file"]
# Received from the harpy wrapper
samplenames = config["samplenames"] 

rule merge_vcfs:
    input: 
        bcf = expand("VariantCall/{sample}.bcf", sample = samplenames),
        index = expand("VariantCall/{sample}.bcf.csi", sample = samplenames)
    output: "VariantCall/variants.raw.bcf"
    log: "VariantCall/variants.raw.stats"
    message: "Merging sample VCFs into single file: {output}"
    default_target: True
    threads: 20
    shell:
        """
        bcftools merge --threads {threads} -o {output} {input.bcf}
        bcftools stats {output} > {log}
        """

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
    output: temp("VariantCall/{sample}.bci")
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
        bc_idx = "VariantCall/{sample}.bci",
        genome = genomefile
    output: vcf = pipe("VariantCall/{sample}.vcf")
    log:  
        runlog = "VariantCall/logs/{sample}.leviathan.log",
        candidates = "VariantCall/logs/{sample}.candidates"
    message: "Calling variants: {wildcards.sample}"
    threads: 3
    shell:
        """
        LEVIATHAN -b {input.bam} -i {input.bc_idx} -g {input.genome} -o {output} -t {threads} --candidates {log.candidates} 2> {log.runlog}
        """

rule vcf2bcf:
    input: "VariantCall/{sample}.vcf"
    output: temp("VariantCall/{sample}.bcf")
    message: "Covnerting to BCF: {input}"
    threads: 1
    shell:        
        """
        bcftools convert -Ob {input} | bcftools sort --output {output}
        """

rule index_bcf:
    input: "VariantCall/{sample}.bcf"
    output: temp("VariantCall/{sample}.bcf.csi")
    message: "Indexing: {input}"
    threads: 1
    shell:
        """
        bcftools index --output {output} {input}
        """