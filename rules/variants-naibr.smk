import os
import re

bam_dir = config["seq_directory"]
samplenames = config["samplenames"] 
extra = config.get("extra", "") 

outdir = "Variants/naibr"

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
        bam = bam_dir + "/{sample}.bam",
        configfile = outdir + "configs/{sample}.config"
    output:
        bedpe     = outdir + "{sample}/{sample}.bedpe",
        bedpe_fmt = outdir + "{sample}/{sample}.reformat.bedpe" 
    threads:
        8        
    params:
        outdir + "{wildcards.sample}"
    message:
        "Calling variants: {wildcards.sample}"
    shell:
        """
        naibr {input.configfile}
        mv {params}/NAIBR.bedpe {output.bedpe}
        mv {params}/NAIBR.reformat.bedpe {output.bedpe_fmt}
        """