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
        "min_mapq" : 30,
        "d"        : 10000,
        "min_sv"   : 1000,
        "k"        : 3
    }
    if args != "":
        words = [i for i in re.split("\s|=", args) if len(i) > 0]
        for i in zip(words[::2], words[1::2]):
            argsDict[i[0]] = i[1]
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
        outdir + "/configs/{sample}.config"
    message:
        "Creating naibr config file: {wildcards.sample}"
    params:
        lambda wc: wc.get("sample")
    run:
        argdict = process_args(extra)
        with open(output[0], "w") as conf:
            _ = conf.write(f"bam_file={input[0]}\n")
            _ = conf.write(f"prefix={params[0]}\n")
            _ = conf.write(f"outdir=Variants/naibr/{params[0]}\n")
            for i in argdict:
                _ = conf.write(f"{i}={argdict[i]}\n")

rule call_sv:
    input:
        bam   = bam_dir + "/{sample}.bam",
        bai   = bam_dir + "/{sample}.bam.bai",
        conf  = outdir + "/configs/{sample}.config"
    output:
        bedpe = outdir + "/{sample}.bedpe",
        refmt = outdir + "/IGV/{sample}.reformat.bedpe",
        fail  = outdir + "/filtered/{sample}.fail.bedpe",
        vcf   = outdir + "/vcf/{sample}.vcf" 
    log:
        outdir + "/logs/{sample}.log"
    threads:
        8        
    params:
        outdir = lambda wc: outdir + "/" + wc.get("sample"),
        sample = lambda wc: wc.get("sample")
    message:
        "Calling variants: {wildcards.sample}"
    log:
        outdir + "log/{sample}.log" 
    shell:
        """
        echo "threads={threads}" >> {input.conf}
        naibr {input.conf} > {log}.tmp 2>&1
        grep -v "pairs/s" {log}.tmp > {log} && rm {log}.tmp
        inferSV.py {params.outdir}/{params.sample}.bedpe -f {output.fail} > {output.bedpe}
        mv {params.outdir}/{params.sample}.reformat.bedpe {output.refmt}
        mv {params.outdir}/{params.sample}.vcf {output.vcf}
        #mv Variants/naibrlog/{params.sample}.log {log}
        rm -rf {params.outdir}
        """

rule link_genome:
    input:
        genomefile
    output: 
        f"Genome/{bn}"
    message:
        "Symlinking {input} to Genome/"
    shell: 
        "ln -sr {input} {output}"

rule faidx_genome:
    input: 
        f"Genome/{bn}"
    output: 
        f"Genome/{bn}.fai"
    message:
        "Indexing {input}"
    log:
        f"Genome/{bn}.faidx.log"
    shell: 
        """
        samtools faidx --fai-idx {output} {input} 2> {log}
        """

rule report:
    input:
        bedpe = outdir + "/{sample}.bedpe",
        fai   = f"Genome/{bn}.fai"
    output:
        outdir + "/reports/{sample}.naibr.html"
    message:
        "Creating report: {wildcards.sample}"
    script:
        "reportNaibr.Rmd"

rule log_runtime:
    output:
        outdir + "/logs/harpy.variants.log"
    message:
        "Creating record of relevant runtime parameters: {output}"
    run:
        argdict = process_args(extra)
        with open(output[0], "w") as f:
            _ = f.write("The harpy variants sv module ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write(f"The directory with alignments: {bam_dir}\n\n")
            _ = f.write("naibr variant calling ran using these configurations:\n")
            _ = f.write(f"    bam_file=BAMFILE\n")
            _ = f.write(f"    prefix=PREFIX\n")
            _ = f.write(f"    outdir=Variants/naibr/PREFIX\n")
            for i in argdict:
                _ = f.write(f"    {i}={argdict[i]}\n")

rule all:
    input:
        expand(outdir + "/{sample}.bedpe",      sample = samplenames),
        expand(outdir + "/reports/{sample}.naibr.html", sample = samplenames),
        outdir + "/logs/harpy.variants.log"
    default_target: True
    message:
        "Variant calling completed!"
    shell:
        "rm -rf Variants/naibrlog"
