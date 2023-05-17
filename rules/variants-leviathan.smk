import os

bam_dir = config["seq_directory"]
genomefile = config["genomefile"]
samplenames = config["samplenames"] 
extra = config.get("extra", "") 

bn = os.path.basename(genomefile)

rule index_alignment:
    input:
        bam_dir + "/{sample}.bam"
    output:
        bam_dir + "/{sample}.bam.bai"
    message:
        "Indexing alignment: {wildcards.sample}"
    benchmark:
        "Benchmark/Variants/leviathan/indexbam.{sample}.txt"
    shell:
        "sambamba index {input} {output} 2> /dev/null"

rule index_barcode:
    input: 
        bam = bam_dir + "/{sample}.bam",
        bai = bam_dir + "/{sample}.bam.bai"
    output:
        temp("Variants/leviathan/lrezIndexed/{sample}.bci")
    message:
        "Indexing barcodes: {wildcards.sample}"
    benchmark:
        "Benchmark/Variants/leviathan/indexbc.{sample}.txt"
    threads: 4
    shell:
        "LRez index bam --threads {threads} -p -b {input.bam} -o {output}"

rule link_genome:
	input:
		genomefile
	output: 
		f"Assembly/{bn}"
	message:
		"Symlinking {input} to Assembly/"
	shell: 
		"ln -sr {input} {output}"

rule index_faidx_genome:
    input: 
        f"Assembly/{bn}"
    output: 
        f"Assembly/{bn}.fai"
    message:
        "Indexing {input}"
    log:
        f"Assembly/{bn}.faidx.log"
    shell: 
        """
        samtools faidx --fai-idx {output} {input} 2> {log}
        """

rule index_bwa_genome:
    input: 
        f"Assembly/{bn}"
    output: 
        multiext(f"Assembly/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    message:
        "Indexing {input}"
    log:
        f"Assembly/{bn}.idx.log"
    shell: 
        """
        bwa index {input} 2> {log}
        """

rule leviathan_variantcall:
    input:
        bam = bam_dir + "/{sample}.bam",
        bai = bam_dir + "/{sample}.bam.bai",
        bc_idx = "Variants/leviathan/lrezIndexed/{sample}.bci",
        genome = f"Assembly/{bn}"
    output:
        pipe("Variants/leviathan/{sample}.vcf")
    log:  
        runlog = "Variants/leviathan/logs/{sample}.leviathan.log",
        candidates = "Variants/leviathan/logs/{sample}.candidates"
    message:
        "Calling variants: {wildcards.sample}"
    benchmark:
        "Benchmark/Variants/leviathan/variantcall.{sample}.txt"
    params:
        extra = extra
    threads: 3
    shell:
        "LEVIATHAN -b {input.bam} -i {input.bc_idx} {params} -g {input.genome} -o {output} -t {threads} --candidates {log.candidates} 2> {log.runlog}"

rule sort_bcf:
    input:
        "Variants/leviathan/{sample}.vcf"
    output:
        "Variants/leviathan/{sample}.bcf"
    message:
        "Sorting and converting to BCF: {wildcards.sample}"
    threads: 1
    params:
        "{wildcards.sample}"
    benchmark:
        "Benchmark/Variants/leviathan/sortbcf.{sample}.txt"
    shell:        
        "bcftools sort -Ob --output {output} {input} 2> /dev/null"

rule sv_stats:
    input: 
        "Variants/leviathan/{sample}.bcf"
    output: 
        "Variants/leviathan/reports/stats/{sample}.sv.stats"
    message:
        "Getting SV stats for {wildcards.sample}"
    benchmark:
        "Benchmark/Variants/leviathan/stats.{sample}.txt"
    threads: 1
    shell:
        """
        echo -e "sample\\tcontig\\tposition_start\\tposition_end\\tlength\\ttype\\tn_barcodes\\tn_pairs" > {output}
        bcftools query -f '{wildcards.sample}\\t%CHROM\\t%POS\\t%END\\t%SVLEN\\t%SVTYPE\\t%BARCODES\\t%PAIRS\\n' {input} >> {output}
        """

rule sv_report:
    input:	
        bcf = "Variants/leviathan/{sample}.bcf",
        statsfile = "Variants/leviathan/reports/stats/{sample}.sv.stats"
    output:	
        "Variants/leviathan/reports/{sample}.SV.html"
    message:
        "Generating SV report: {wildcards.sample}"
    script:
        "reportLeviathan.Rmd"

rule all_bcfs:
    input: 
        bcf = expand("Variants/leviathan/{sample}.bcf", sample = samplenames),
        reports = expand("Variants/leviathan/reports/{sample}.SV.html", sample = samplenames)
    default_target: True
    message:
        "Variant calling is complete!"