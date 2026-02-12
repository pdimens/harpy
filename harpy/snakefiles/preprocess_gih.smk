import os
import re
import shutil
from harpy.common.preprocess import needs_stagger

wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+"

VERSION=4.0
fqlist       = config["Inputs"]
skip_reports = config["Workflow"]["reports"]["skip"]
#bc_len       = config["Parameters"]["barcode-length"]
insert_min   = config["Parameters"]["minimum-length"]
me_seq       = config["Parameters"]["ME-sequence"] 
overlap      = config["Parameters"]["ME-overlap"] 
has_pigz     = bool(shutil.which('pigz'))

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
    input: collect("stagger/{sample}.R1.fq.gz", sample = samplenames)


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

rule barcode_padding:
    input:
        FQ1 = get_fq1,
        FQ2 = get_fq2,
        info = "ME_position/{sample}.info",
        summary = "ME_position/{sample}.info.summ"
    output:
        FQ1 = "stagger/{sample}.R1.fq.gz",
        FQ2 = "stagger/{sample}.R2.fq.gz"
    log:
        "logs/{sample}.stagger.log"
    threads:
        4 if has_pigz else 1
    run:
        if needs_stagger(input.summary):
            shell(f"stagger-GIH -t {threads} -b 20000 stagger/{wildcards.sample} {input.info} {input.FQ1} {input.FQ2} 2> {log[0]}")
        else:
            shell(f"ln -sr {input.FQ1} {output.FQ1} 2> {log[0]}")
            shell(f"ln -sr {input.FQ2} {output.FQ2} 2>> {log[0]}")

rule extract_barcodes:
    input:
        FQ1 = "stagger/{sample}.R1.fq.gz",
        FQ2 = "stagger/{sample}.R1.fq.gz",
        pheniqs_conf = "workflow/pheniqs_config.json"
    output:
        FW = "extract/{sample}.R1.bam",
        RV = "extract/{sample}.R2.bam"
    log:
        "logs/extract/{sample}.json"
    shell:
        """
        pheniqs mux \
            --input {input.FQ1} \
            --input {input.FQ2} \
            --output {output.FW} \
            --output {output.RV} \
            -c {input.pheniqs_conf} \
            --quality \
            --report {log}
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