import os
import re
import logging

onstart:
    logfile_handler = logger_manager._default_filehandler(config["Workflow"]["snakemake"]["log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

fqlist        = config["Inputs"]
skip_reports  = config["Workflow"]["reports"]["skip"]
bc_len        = config["Parameters"]["barcode_length"]
insert_min    = config["Parameters"]["minimum_length"]

bn_r = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"

def get_fq1(wildcards):
    # returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    r = re.compile(fr"(.*/{re.escape(wildcards.sample)})([_\.]1|[_\.]F|[_\.]R1(?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$", flags = re.IGNORECASE)
    return list(filter(r.match, fqlist))

def get_fq2(wildcards):
    # returns a list of fastq files for read 2 based on *wildcards.sample*, e.g.
    r = re.compile(fr"(.*/{re.escape(wildcards.sample)})([_\.]2|[_\.]R|[_\.]R2(?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$", flags = re.IGNORECASE)
    return list(filter(r.match, fqlist))

rule qc:
    input:
        fw   = get_fq1,
        rv   = get_fq2
    output:
        fw = ,
        fw_short = ,
        rv = ,
        rv_short =
    log:
        "logs/qc/{sample}_precheck-inserts.log"
    params:
        "-a CTGTCTCTTATACACATCT",
        "-A AGATCGGAAGAGC",
        "--nextseq-trim=20",
        "--match-read-wildcards",
        "--action trim",
        f"--minimum-length {bc_len + insert_min}"
    threads:
        10
    conda:

    container:

    shell:
        """
        cutadapt {params} \
            --too-short-output {output.fw_short} \
            --too-short-paired-output {output.rv_short} \
            --cores {threads} \
            --output {output.fw} \
            --paired-output {output.rv} \
            {input} > {log} 
        """

rule findMEseq:

    params:
        "AGATGTGTATAAGAGACAG"


rule determine_stagger:
    input:
        FQ1 = ,
        FQ2 = ,
        summary = 
    shell:
        """
        include_stagger=\$(python $baseDir/stagger_check.py {input.summary})

        # Determine bead complexity
        ln -s $baseDir/ .
        samtools view {input.FQ1} | awk 'NR>100000 && NR<=2000000 {n=split(\$0,arr,"CB:Z:"); if(n>1){split(arr[2],result,"\t"); print result[1];}}' > temp_barcodes.txt
        complexity=\$(python $baseDir/determine_bead_complexity.py --include_stagger "\$include_stagger" --temp_barcodes temp_barcodes.txt)
        grep "\$complexity" 12nt-barcodes-with-stagger-and-bead-complexity.txt > 12nt-barcodes-subset.txt

        # Convert 43nt barcode sequences to ‘haplotagging’ format AxxCxxBxxDxx for both R1 and R2 files.
        awk -v stagger="\$include_stagger" -f $baseDir/convert_bam_to_ACBD_format.awk 12nt-barcodes-subset.txt <(samtools view {input.FQ1}) | pigz -p 8 > \${sample}_converted.R1.fastq.gz
        mv barcode_log.log \${sample}_converted_R1_barcode_log.log
        awk -v stagger="\$include_stagger" -f $baseDir/convert_bam_to_ACBD_format.awk 12nt-barcodes-subset.txt <(samtools view {input.FQ1}) | pigz -p 8 > \${sample}_converted.R2.fastq.gz
        mv barcode_log.log \${sample}_converted_R2_barcode_log.log
        """

rule pheniqs_demux:
    input:
        FQ1 = ,
        FQ2 = ,
        pheniqs_conf = "workflow/pheniqs_config.json"
    output:
        FQ1 = "pheniqs/{sample}_extract_R1.bam",
        FQ2 = "pheniqs/{sample}_extract_R2.bam"
    log:
        "logs/pheniqs/{sample}.json"
    shell:
        """
        pheniqs mux \
            --input {input.FQ1} \
            --input {input.FQ2} \
            --output {output.FQ1} \
            --output {output.FQ2} \
            -c {input.pheniqs_conf} \
            --quality \
            --report {log}
        """

