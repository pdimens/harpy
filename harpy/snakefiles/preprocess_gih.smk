import os
import re
import shutil
from harpy.common.preprocess import needs_stagger, padUMI

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

rule find_ME_sequence:
    input:
        get_fq1
    output:
        tmp = "ME_position/{sample}.info",
        summ = temp("ME_position/{sample}.info.summ")
    log:
        "logs/findME/{sample}.log"
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
            cutadapt {params.static} --info-file {output.tmp} --cores {threads} {input}
            {params.awk} {output.tmp} > {output.summ}
        }} 2> {log}
        """

rule all:
    default_target: True
    #input: collect("ME_position/{sample}.info", sample = samplenames)
    input: collect("stagger/{sample}.R1.fq.gz", sample = samplenames)


rule pad_UMI:
    input:
        fq = get_fq1,
        info = "ME_position/{sample}.info",
        summary = "ME_position/{sample}.info.summ"
    output:
        fq = "stagger/{sample}.R1.fq.gz"
    log:
        "logs/{sample}.stagger.log"
    threads:
        4 if has_pigz else 1
    run:
        try:
            result = needs_stagger(input.summary)
            if result:
                padUMI(input.info, input.fq, output.fq, log[0], threads, 20000)
            else:
                shell(f"ln -sr {input.fq} {output.fq}")
        except Exception as e:
            with open(logs[0], "w") as _log:
                _log.write(f"{e}")

#rule move_barcodes:
#    input:
#        FQ1 = ,
#        FQ2 = ,
#        pheniqs_conf = "workflow/pheniqs_config.json"
#    output:
#        FQ1 = "pheniqs/{sample}_extract_R1.bam",
#        FQ2 = "pheniqs/{sample}_extract_R2.bam"
#    log:
#        "logs/pheniqs/{sample}.json"
#    shell:
#        """
#        pheniqs mux \
#            --input {input.FQ1} \
#            --input {input.FQ2} \
#            --output {output.FQ1} \
#            --output {output.FQ2} \
#            -c {input.pheniqs_conf} \
#            --quality \
#            --report {log}
#        """
#
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