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
    output: vcf = temp("Variants/leviathan/{sample}.vcf")
    log:  
        runlog = "Variants/leviathan/logs/{sample}.leviathan.log",
        candidates = "Variants/leviathan/logs/{sample}.candidates"
    message: "Calling variants: {wildcards.sample}"
    params:
        extra = extra
    threads: 4
    shell:
        """
        LEVIATHAN -b {input.bam} -i {input.bc_idx} {params} -g {input.genome} -o {output} -t {threads} --candidates {log.candidates} 2> {log.runlog}
        """

rule add_samplename:
    input: "Variants/leviathan/{sample}.vcf"
    output: 
        bcf = temp("Variants/leviathan/{sample}.bcf"),
        namefile = temp(".{sample}.name")
    message: "Covnerting to BCF: {input}"
    threads: 1
    params: "{wildcards.sample}"
    shell:        
        """
        echo {params} > {output.namefile}
        bcftools reheader --samples {output.namefile} {input} | bcftools sort -Ob --output {output.bcf} 2> /dev/null
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