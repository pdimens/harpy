import sys
import os
import re

bam_dir     = config["seq_directory"]
samplenames = config["samplenames"] 
extra       = config.get("extra", "") 
groupfile   = config["groupings"]
genomefile  = config["genomefile"]
bn          = os.path.basename(genomefile)
outdir      = "Variants/naibr-pop"

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

# create dictionary of population => filenames
## this makes it easier to set the snakemake rules/wildcards
## exits with an error if the groupfile has samples not in the bam folder
def pop_manifest(infile, dirn, sampnames):
	d = dict()
	absent = []
	with open(infile) as f:
		for line in f:
			samp, pop = line.rstrip().split()
			if samp not in sampnames:
				absent.append(samp)
			samp = f"{dirn}/{samp}.bam"
			if pop not in d.keys():
				d[pop] = [samp]
			else:
				d[pop].append(samp)
	if absent:
		sys.tracebacklimit = 0
		raise ValueError(f"{len(absent)} sample(s) in \033[1m{infile}\033[0m not found in \033[1m{dirn}\033[0m directory:\n\033[33m" + ", ".join(absent) + "\033[0m")
	return d

popdict     = pop_manifest(groupfile, bam_dir, samplenames)
populations = popdict.keys()

rule bamlist:
	output:
		expand(outdir + "/input/{pop}.list", pop = populations)
	message:
		"Creating file lists for each population."
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
		temp(outdir + "/input/{population}.bam")
	message:
		"Merging alignments: Population {wildcards.population}"
	shell:
		"samtools merge -b {input} -o {output}"

rule index_merged:
	input:
		outdir + "/input/{population}.bam"
	output:
		temp(outdir + "/input/{population}.bam.bai")
	message:
		"Indexing merged alignments: Population {wildcards.population}"
	wildcard_constraints:
		population = "[a-zA-Z0-9_-]*"
	shell:
		"sambamba index {input} {output} 2> /dev/null"

rule create_config:
    input:
        outdir + "/input/{population}.bam"
    output:
        temp(outdir + "/configs/{population}.config")
    message:
        "Creating naibr config file: {wildcards.population}"
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
        bam    = outdir + "/input/{population}.bam",
        bai    = outdir + "/input/{population}.bam.bai",
        config = outdir + "/configs/{population}.config"
    output:
        bedpe  = outdir + "/{population}.bedpe",
        refmt  = outdir + "/IGV/{population}.reformat.bedpe",
        fail   = outdir + "/filtered/{population}.fail.bedpe",
        vcf    = outdir + "/vcf/{population}.vcf"
    threads:
        8        
    params:
        outdir + "{wildcards.population}"
    message:
        "Calling variants: {wildcards.population}"
    log:
        outdir + "logs/{population}.log",
    shell:
        """
        naibr {input.configfile} 2>&1 > {log}
        inferSV.py {params}/NAIBR.bedpe -f {output.fail} > {output.bedpe}
        mv {params}/NAIBR.reformat.bedpe {output.refmt}
        mv {params}/NAIBR.vcf {output.vcf}
        rm -rf {params}
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

rule report:
    input:
        fai   = f"Assembly/{bn}.fai",
        bedpe = outdir + "/{population}.bedpe"
    output:
        outdir + "/reports/{population}.naibr.html"
    message:
        "Creating report: {wildcards.population}"
	script:
		"reportNaibr.Rmd"

rule report_pop:
    input:
        fai   = f"Assembly/{bn}.fai",
        bedpe = expand(outdir + "/{pop}.bedpe", pop = populations)
    output:
        outdir + "/reports/naibr.summary.html"
    message:
        "Creating report: {wildcards.population}"
	script:
		"reportNaibrPop.Rmd"

rule all:
    input:
        expand(outdir + "/{pop}.bedpe",      pop = populations),
        expand(outdir + "/{pop}.naibr.html", pop = populations).
        outdir + "/reports/naibr.summary.html"
    default_target: True
    message:
        "Variant calling completed!"