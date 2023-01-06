#BUG leviathan doesn't seem to output sample ID into the vcf
bam_dir = config["seq_directory"]
genomefile = config["genomefile"]
samplenames = config["samplenames"] 

rule merge_bcfs:
    input: 
        bcf = expand("Variants/leviathan/{sample}.bcf", sample = samplenames),
        index = expand("Variants/leviathan/{sample}.bcf.csi", sample = samplenames)
    output: "Variants/leviathan/variants.raw.bcf"
    log: "Variants/leviathan/variants.raw.stats"
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
    output: temp("Variants/leviathan/{sample}.bci")
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
        bc_idx = "Variants/leviathan/{sample}.bci",
        genome = genomefile
    output: vcf = pipe("Variants/leviathan/{sample}.vcf")
    log:  
        runlog = "Variants/leviathan/logs/{sample}.leviathan.log",
        candidates = "Variants/leviathan/logs/{sample}.candidates"
    message: "Calling variants: {wildcards.sample}"
    threads: 3
    shell:
        """
        LEVIATHAN -b {input.bam} -i {input.bc_idx} -g {input.genome} -o {output} -t {threads} --candidates {log.candidates} 2> {log.runlog}
        """

rule vcf2bcf:
    input: "Variants/leviathan/{sample}.vcf"
    output: temp("Variants/leviathan/{sample}.bcf")
    message: "Covnerting to BCF: {input}"
    threads: 1
    shell:        
        """
        bcftools convert -Ob {input} | bcftools sort --output {output}
        """

rule index_bcf:
    input: "Variants/leviathan/{sample}.bcf"
    output: temp("Variants/leviathan/{sample}.bcf.csi")
    message: "Indexing: {input}"
    threads: 1
    shell:
        """
        bcftools index --output {output} {input}
        """