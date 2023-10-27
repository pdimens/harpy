import sys
import os
import re

bam_dir     = config["seq_directory"]
samplenames = config["samplenames"] 
extra       = config.get("extra", "") 
groupfile   = config["groupings"]
genomefile  = config["genomefile"]
molecule_distance = config["molecule_distance"]
bn          = os.path.basename(genomefile)
outdir      = "Variants/naibr-pop"
genome_zip  = True if bn.lower().endswith(".gz") else False
if genome_zip:
    bn = bn[:-3]

def process_args(args):
    argsDict = {
        "min_mapq" : 30,
        "d"        : molecule_distance,
        "min_sv"   : 1000,
        "k"        : 3,
    }
    if args != "":
        words = [i for i in re.split("\s|=", args) if len(i) > 0]
        for i in zip(words[::2], words[1::2]):
            argsDict[i[0]] = i[1]
    return argsDict

# create dictionary of population => filenames
## this makes it easier to set the snakemake rules/wildcards
## exits with an error if the groupfile has samples not in the bam folder
def pop_manifest(infile, dirn, sampnames):
    d = dict()
    absent = []
    with open(infile) as f:
        for line in f:
            samp, pop = line.rstrip().split()
            if samp.lstrip().startswith("#"):
                continue
            if samp not in sampnames:
                absent.append(samp)
            samp = f"{dirn}/{samp}.bam"
            if pop not in d.keys():
                d[pop] = [samp]
            else:
                d[pop].append(samp)
    if absent:
        sys.tracebacklimit = 0
        raise ValueError(f"{len(absent)} sample(s) in \033[1m{infile}\033[0m not found in \033[1m{dirn}\033[0m directory:\n\033[33m" + ", ".join(absent) + "\033[0m" + "\n")
    return d

popdict     = pop_manifest(groupfile, bam_dir, samplenames)
populations = popdict.keys()

rule copy_groupings:
    input:
        groupfile
    output:
        outdir + "/logs/sample.groups"
    message:
        "Logging {input}"
    run:
        with open(input[0], "r") as infile, open(output[0], "w") as outfile:
            _ = [outfile.write(i) for i in infile.readlines() if not i.lstrip().startswith("#")]

rule bamlist:
    input:
        outdir + "/logs/sample.groups"
    output:
        expand(outdir + "/input/{pop}.list", pop = populations)
    message:
        "Creating population file lists."
    run:
        for p in populations:
            bamlist = popdict[p]
            with open(f"{outdir}/input/{p}.list", "w") as fout:
                for bamfile in bamlist:
                    _ = fout.write(bamfile + "\n")

rule merge_populations:
    input: 
        bamlist  = outdir + "/input/{population}.list",
        bamfiles = lambda wc: expand("{sample}", sample = popdict[wc.population]) 
    output:
        bam = temp(outdir + "/input/{population}.bam"),
        bai = temp(outdir + "/input/{population}.bam.bai")
    message:
        "Merging alignments: Population {wildcards.population}"
    threads:
        2
    shell:
        "samtools merge -o {output.bam}##idx##{output.bai} --threads {threads} --write-index -b {input.bamlist}"

rule create_config:
    input:
        outdir + "/input/{population}.bam"
    output:
        outdir + "/configs/{population}.config"
    message:
        "Creating naibr config file: {wildcards.population}"
    params:
        lambda wc: wc.get("population")
    run:
        argdict = process_args(extra)
        with open(output[0], "w") as conf:
            _ = conf.write(f"bam_file={input[0]}\n")
            _ = conf.write(f"outdir=Variants/naibr-pop/{params[0]}\n")
            _ = conf.write(f"prefix={params[0]}\n")
            for i in argdict:
                _ = conf.write(f"{i}={argdict[i]}\n")

rule call_sv:
    input:
        bam   = outdir + "/input/{population}.bam",
        bai   = outdir + "/input/{population}.bam.bai",
        conf  = outdir + "/configs/{population}.config"
    output:
        bedpe = outdir + "/{population}.bedpe",
        refmt = outdir + "/IGV/{population}.reformat.bedpe",
        fail  = outdir + "/bad_candidates/{population}.fail.bedpe",
        vcf   = outdir + "/vcf/{population}.vcf"
    log:
        outdir + "/logs/{population}.log"
    threads:
        8        
    params:
        population = lambda wc: wc.get("population"),
        outdir     = lambda wc: outdir + "/" + wc.get("population")
    message:
        "Calling variants: {wildcards.population}"
    shell:
        """
        if ! grep -q "threads" {input.conf}; then
            echo "threads={threads}" >> {input.conf}
        fi
        naibr {input.conf} > {log}.tmp 2>&1
        grep -v "pairs/s" {log}.tmp > {log} && rm {log}.tmp
        inferSV.py {params.outdir}/{params.population}.bedpe -f {output.fail} > {output.bedpe}
        mv {params.outdir}/{params.population}.reformat.bedpe {output.refmt}
        mv {params.outdir}/{params.population}.vcf {output.vcf}
        #mv Variants/naibrlog/{params.population}.log {log}
        rm -rf {params.outdir}
        """

rule genome_link:
    input:
        genomefile
    output: 
        f"Genome/{bn}"
    message: 
        "Symlinking {input}"
    shell: 
        """
        if (file {input} | grep -q compressed ) ;then
            # is regular gzipped, needs to be BGzipped
            zcat {input} | bgzip -c > {output}
        elif (file {input} | grep -q BGZF ); then
            # is bgzipped, just linked
            ln -sr {input} {output}
        else
            # isn't compressed, just linked
            ln -sr {input} {output}
        fi
        """

if genome_zip:
    rule genome_compressed_faidx:
        input: 
            f"Genome/{bn}"
        output: 
            gzi = f"Genome/{bn}.gzi",
            fai = f"Genome/{bn}.fai"
        message:
            "Indexing {input}"
        log:
            f"Genome/{bn}.faidx.gzi.log"
        shell: 
            "samtools faidx --gzi-idx {output.gzi} --fai-idx {output.fai} {input} 2> {log}"
else:
    rule genome_faidx:
        input: 
            f"Genome/{bn}"
        output: 
            f"Genome/{bn}.fai"
        message:
            "Indexing {input}"
        log:
            f"Genome/{bn}.faidx.log"
        shell:
            "samtools faidx --fai-idx {output} {input} 2> {log}"

rule report:
    input:
        fai   = f"Genome/{bn}.fai",
        bedpe = outdir + "/{population}.bedpe"
    output:
        outdir + "/reports/{population}.naibr.html"
    message:
        "Creating report: {wildcards.population}"
    script:
        "reportNaibr.Rmd"

rule report_pop:
    input:
        fai   = f"Genome/{bn}.fai",
        bedpe = expand(outdir + "/{pop}.bedpe", pop = populations)
    output:
        outdir + "/reports/naibr.pop.summary.html"
    message:
        "Creating summary report"
    script:
        "reportNaibrPop.Rmd"

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
            _ = f.write(f"The directory with alignments: {bam_dir}\n")
            _ = f.write(f"The sample grouping file: {groupfile}\n\n")
            _ = f.write("naibr variant calling ran using these configurations:\n")
            _ = f.write(f"    bam_file=BAMFILE\n")
            _ = f.write(f"    prefix=PREFIX\n")
            _ = f.write(f"    outdir=Variants/naibr/PREFIX\n")
            for i in argdict:
                _ = f.write(f"    {i}={argdict[i]}\n")

rule all:
    default_target: True
    input:
        expand(outdir + "/{pop}.bedpe",      pop = populations),
        expand(outdir + "/reports/{pop}.naibr.html", pop = populations),
        outdir + "/reports/naibr.pop.summary.html",
        outdir + "/logs/harpy.variants.log"
    message:
        "Variant calling completed!"
    shell:
        "rm -rf Variants/naibrlog"