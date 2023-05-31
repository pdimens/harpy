import os
import re

bam_dir     = config["seq_directory"]
samplenames = config["samplenames"] 
extra       = config.get("extra", "") 
genomefile  = config["genomefile"]
bn          = os.path.basename(genomefile)
outdir      = "Variants/naibr"

def process_args(args):
    argsDict = {
        min_mapq : 30
        d        : 10000
        min_sv   : 1000
        k        : 3
    }
    if args != "":
        words = [i for i in re.split("\s|=", args) if len(i) > 0]
        for i in zip(words[::2], words[1::2]):
            argsDict[i[0]]] = i[1]
    return argsDict

rule index_alignment:
    input:
        bam_dir + "/{sample}.bam"
    output:
        bam_dir + "/{sample}.bam.bai"
    message:
        "Indexing alignment: {wildcards.sample}"
    shell:
        "sambamba index {input} {output} 2> /dev/null"

rule create_config:
    input:
        bam_dir + "/{sample}.bam"
    output:
        temp(outdir + "configs/{sample}.config")
    message:
        "Creating naibr config file: {wildcards.sample}"
    params:
        extra
    run:
        from multiprocessing import cpu_count
        argdict = process_args(params)
        with open(output[0], "w") as conf:
            _ = conf.write(f"bam_file={wildcards.sample}\n")
            _ = conf.write(f"outdir={wildcards.sample}\n")
            _ = conf.write(f"threads={cpu_count()}\n")
            for i in argdict:
                _ = conf.write(f"{i}={argdict[i]}\n")

rule call_sv:
    input:
        bam        = bam_dir + "/{sample}.bam",
        configfile = outdir + "configs/{sample}.config"
    output:
        bedpe      = outdir + "{sample}/{sample}.bedpe",
        bedpe_fmt  = outdir + "{sample}/{sample}.reformat.bedpe" 
        bedpe_fail = outdir + "{sample}/{sample}.fail.bedpe"
        vcf        = outdir + "{sample}/{sample}.vcf" 
    threads:
        8        
    params:
        outdir + "{wildcards.sample}"
    message:
        "Calling variants: {wildcards.sample}"
    log:
        outdir + "{sample}/{sample}.log" 
    shell:
        """
        naibr {input.configfile} 2>&1 > {log}
        inferSV.py {params}/NAIBR.bedpe -f > {output.bedpe}
        mv {params}/NAIBR.reformat.bedpe {output.bedpe_fmt}
        mv {params}/NAIBR.vcf {output.vcf}
        """

rule link_genome:
	input:
		genomefile
	output: 
		f"Assembly/{bn}"
	message:
		"Symlinking {input} to Assembly/"
	shell: 
		"ln -sr {input} {output}"

rule faidx_genome:
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

rule report:
    input:
        bedpe = outdir + "{sample}/{sample}.bedpe",
        fai   = f"Assembly/{bn}.fai"
    output:
        outdir + "{sample}/{sample}.naibr.html"
    message:
        "Creating report: {wildcards.sample}"
	script:
		"reportNaibr.Rmd"

rule all:
    input:
        expand(outdir + "{sample}/{sample}.bedpe", sample = samplenames),
        expand(outdir + "{sample}/{sample}.naibr.html", sample = samplenames)
    default_target: True
    message:
        "Variant calling completed!"
