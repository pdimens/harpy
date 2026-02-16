import os
import re

wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

VERSION=4.0
fqlist       = config["Inputs"]
skip_reports = config["Workflow"]["reports"]["skip"]
#bc_len       = config["Parameters"]["barcode-length"]
insert_min   = config["Parameters"]["minimum-length"]
me_seq       = config["Parameters"]["ME-sequence"] 
overlap      = config["Parameters"]["ME-overlap"] 

bn_r = r"([_\.][12]|[_\.][FR]|[_\.]R[12](?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$"
samplenames = {re.sub(bn_r, "", os.path.basename(i), flags = re.IGNORECASE) for i in fqlist}

def get_fq1(wildcards):
    # returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    r = re.compile(fr"(.*/{re.escape(wildcards.sample)})([_\.]1|[_\.]F|[_\.]R1(?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$", flags = re.IGNORECASE)
    return list(filter(r.match, fqlist))

def get_fq2(wildcards):
    # returns a list of fastq files for read 2 based on *wildcards.sample*, e.g.
    r = re.compile(fr"(.*/{re.escape(wildcards.sample)})([_\.]2|[_\.]R|[_\.]R2(?:\_00[0-9])*)?\.((fastq|fq)(\.gz)?)$", flags = re.IGNORECASE)
    return list(filter(r.match, fqlist))


rule all:
    default_target: True
    input: collect("stagger/{sample}.stagger.bam", sample = samplenames)


rule find_ME_seq:
    input:
        get_fq1
    output:
        info = temp("ME_position/{sample}.info"),
        summ = temp("ME_position/{sample}.info.summ")
    log:
        "logs/find_ME_seq/{sample}.log"
    params:
        static = f"-g {me_seq} --overlap {overlap} -e 0.11 --match-read-wildcards --action none -o /dev/null",
        awk = """awk -F '\\t' '{a[$2]++; if($2>=0) {b[$3]++; c[$4]++;} else next;} END {print "col2=mismatch"; for(i in a) print a[i],i; print "\\ncol3=startpost"; for(j in b) print b[j],j; print "\\ncol4=endpos"; for(k in c) print c[k],k;}'"""
    threads:
        4
    conda:
        "envs/qc.yaml"
    container:
        f"docker://pdimens/harpy:qc_{VERSION}"
    shell:
        """
        {{
            cutadapt {params.static} --info-file {output.info} --cores {threads} {input}
            {params.awk} {output.info} > {output.summ}
        }} 2> {log}
        """

rule pad_barcodes:
    input:
        info = "ME_position/{sample}.info",
        summary = "ME_position/{sample}.info.summ",
        FQ1 = get_fq1,
        FQ2 = get_fq2
    output:
        temp("stagger/{sample}.stagger.bam")
    log:
        "logs/{sample}.stagger.log"
    threads:
        2
    shell:
        "stagger-GIH {input} | samtools import -s - > {output} 2> {log}"

rule extract_barcodes:
    input:
        stagger = "stagger/{sample}.stagger.bam",
        pheniqs_conf = "workflow/pheniqs_config.json"
    output:
        FW = "extract/{sample}.R1.bam",
        RV = "extract/{sample}.R2.bam"
    log:
        json = "logs/extract/{sample}.json",
        err = "logs/{sample}.pheniqs"
    conda:
        "envs/preprocess.yaml"
    container:
        f"docker://pdimens/harpy:preprocess_{VERSION}"
    shell:
        """
        pheniqs mux \
            --input {input.stagger} \
            --input {input.stagger} \
            --output {output.FW} \
            --output {output.RV} \
            -c {input.pheniqs_conf} \
            --quality \
            --report {log.json} 2> {log.err}
        """

#rule convert_to_fastq:
#    input:
#        FQ1 = ,
#        FQ2 = ,
#        summary = 
#    shell:
#        """
#        include_stagger=\$(python $baseDir/stagger_check.py {input.summary})
#
#        # Determine bead complexity
#        ln -s $baseDir/ .
#        samtools view {input.FQ1} | awk 'NR>100000 && NR<=2000000 {n=split(\$0,arr,"CB:Z:"); if(n>1){split(arr[2],result,"\t"); print result[1];}}' > temp_barcodes.txt
#        complexity=\$(python $baseDir/determine_bead_complexity.py --include_stagger "\$include_stagger" --temp_barcodes temp_barcodes.txt)
#        grep "\$complexity" 12nt-barcodes-with-stagger-and-bead-complexity.txt > 12nt-barcodes-subset.txt
#
#        # Convert 43nt barcode sequences to ‘haplotagging’ format AxxCxxBxxDxx for both R1 and R2 files.
#        awk -v stagger="\$include_stagger" -f $baseDir/convert_bam_to_ACBD_format.awk 12nt-barcodes-subset.txt <(samtools view {input.FQ1}) | pigz -p 8 > \${sample}_converted.R1.fastq.gz
#        mv barcode_log.log \${sample}_converted_R1_barcode_log.log
#        awk -v stagger="\$include_stagger" -f $baseDir/convert_bam_to_ACBD_format.awk 12nt-barcodes-subset.txt <(samtools view {input.FQ1}) | pigz -p 8 > \${sample}_converted.R2.fastq.gz
#        mv barcode_log.log \${sample}_converted_R2_barcode_log.log
#        """