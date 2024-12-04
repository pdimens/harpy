import sys
import os
import re

bam_dir     = config["seq_directory"]
samplenames = config["samplenames"] 
extra       = config.get("extra", "") 
groupfile   = config["groupings"]
genomefile  = config["genomefile"]
molecule_distance = config["molecule_distance"]
outdir      = "Variants/naibr-pop"

bn = os.path.basename(genomefile)
if bn.lower().endswith(".gz"):
    validgenome = bn[:-3]
else:
    validgenome = bn

vcffile = config["vcf"]
if vcffile.lower().endswith("bcf"):
    vcfindex = vcffile + ".csi"
else:
    vcfindex = vcffile + ".tbi"

def process_args(args):
    argsDict = {
        "min_mapq" : 30,
        "d"        : molecule_distance,
        "min_sv"   : 1000,
        "k"        : 3,
    }
    if args:
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
            samp = f"{dirn}/phasedbam/{samp}.bam"
            if pop not in d.keys():
                d[pop] = [samp]
            else:
                d[pop].append(samp)
    if absent:
        sys.tracebacklimit = 0
        raise ValueError(f"{len(absent)} sample(s) in \033[1m{infile}\033[0m not found in \033[1m{dirn}\033[0m directory:\n\033[33m" + ", ".join(absent) + "\033[0m" + "\n")
    return d

popdict     = pop_manifest(groupfile, outdir, samplenames)
populations = popdict.keys()

rule genome_link:
    input:
        genomefile
    output: 
        f"Genome/{validgenome}"
    message: 
        "Preprocessing {input}"
    shell: 
        """
        if (file {input} | grep -q compressed ) ;then
            # decompress gzipped
            seqtk seq {input} > {output}
        elif (file {input} | grep -q BGZF ); then
            # decompress bgzipped
            seqtk seq {input} > {output}
        else
            # linked uncompressed
            ln -sr {input} {output}
        fi
        """

rule genome_faidx:
    input: 
        f"Genome/{validgenome}"
    output: 
        f"Genome/{validgenome}.fai"
    message:
        "Indexing {input}"
    log:
        f"Genome/{validgenome}.faidx.log"
    shell:
        "samtools faidx --fai-idx {output} {input} 2> {log}"

rule index_bcf:
    input:
        vcffile
    output:
        vcffile + ".csi"
    message:
        "Indexing {input}"
    shell:
        "bcftools index {input}"

rule index_vcfgz:
    input:
        vcffile
    output:
        vcffile + ".tbi"
    message:
        "Indexing {input}"
    shell:
        "tabix {input}"

rule phase_alignments:
    input:
        vcfindex,
        bam_dir + "/{sample}.bam.bai",
        f"Genome/{validgenome}.fai",
        vcf = vcffile,
        aln = bam_dir + "/{sample}.bam",
        ref = f"Genome/{validgenome}"
    output:
        outdir + "/phasedbam/{sample}.bam"
    message:
        "Phasing: {input.aln}"
    params:
        extra = lambda wc: "--ignore-read-groups --sample " + wc.get("sample") + " --tag-supplementary"
    log:
        outdir + "/logs/whatshap-haplotag/{sample}.phase.log"
    threads:
        4
    shell:
        "whatshap haplotag {params} --output-threads={threads} -o {output} --reference {input.ref} {input.vcf} {input.aln} 2> {log}"

rule log_phasing:
    input:
        expand(outdir + "/logs/whatshap-haplotag/{sample}.phase.log", sample = samplenames)
    output:
        outdir + "/logs/whatshap-haplotag/phasing.log"
    message:
        "Creating log of alignment phasing"
    shell:
        """
        echo -e "sample\\ttotal_alignments\\tphased_alignments" > {output}
        for i in {input}; do
            SAMP=$(basename $i .phaselog)
            echo -e "${{SAMP}}\\t$(grep "Total alignments" $i)\\t$(grep "could be tagged" $i)" |
                sed 's/ \+ /\\t/g' | cut -f1,3,5 >> {output}
        done
        """

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

rule report:
    input:
        fai   = f"Genome/{validgenome}.fai",
        bedpe = outdir + "/{population}.bedpe"
    output:
        outdir + "/reports/{population}.naibr.html"
    message:
        "Creating report: {wildcards.population}"
    script:
        "reportNaibr.Rmd"

rule report_pop:
    input:
        fai   = f"Genome/{validgenome}.fai",
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
            _ = f.write("The alignment files were phased using:\n")
            _ = f.write(f"    whatshap haplotag --reference genome.fasta --linked-read-distance-cutoff {molecule_distance} --ignore-read-groups --tag-supplementary --sample sample_x file.vcf sample_x.bam\n")
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