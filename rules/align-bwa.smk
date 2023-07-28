import os
import re
import glob

seq_dir		= config["seq_directory"]
genomefile 	= config["genomefile"]
extra 		= config.get("extra", "") 
bn 			= os.path.basename(genomefile)
outdir      = "Align/bwa"
## deprecated ##
#Rsep 		= config["Rsep"]
#fqext 		= config["fqext"]
#samplenames = config["samplenames"]

#flist = os.listdir(seq_dir)
flist = [os.path.basename(i) for i in glob.iglob(f"{seq_dir}/*") if not os.path.isdir(i)]
r = re.compile(".*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
fqlist = list(filter(r.match, flist))
bn_r = r"[\.\_][RF](?:[12])?(?:\_00[1-9])*\.f(?:ast)?q(?:\.gz)?$"
samplenames = set([re.sub(bn_r, "", i, flags = re.IGNORECASE) for i in fqlist])
d = dict()
for i in samplenames:
    d[i] = i

def get_fq1(wildcards):
    # code that returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    lst = sorted(glob.glob(seq_dir + "/" + wildcards.sample + "*"))
    r = re.compile(".*[\_\.][FR][1]?(?:\_00[0-9])*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
    fqlist = list(filter(r.match, lst))
    return fqlist

def get_fq2(wildcards):
    # code that returns a list of fastq files for read 2 based on *wildcards.sample*, e.g.
    lst = sorted(glob.glob(seq_dir + "/" + wildcards.sample + "*"))
    r = re.compile(".*[\_\.][R][2]?(?:\_00[0-9])*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
    fqlist = list(filter(r.match, lst))
    return fqlist

rule all:
    input: 
        expand(outdir + "/{sample}.bam", sample = samplenames),
        expand(outdir + "/stats/coverage/{sample}.cov.html", sample = samplenames),
        expand(outdir + "/stats/BXstats/{sample}.bxstats.html", sample = samplenames),
        outdir + "/stats/bwa.stats.html",
        outdir + "/logs/harpy.align.log"
    message:
        "Read mapping completed!"
    default_target: True

rule link_genome:
    input:
        genomefile
    output: 
        f"Assembly/{bn}"
    message: 
        "Symlinking {input}"
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
        "samtools faidx --fai-idx {output} {input} 2> {log}"

rule index_bwa_genome:
    input: 
        f"Assembly/{bn}"
    output: 
        multiext(f"Assembly/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    message:
        "Indexing {input}"
    log:
        f"Assembly/{bn}.idx.log"
    shell: 
        "bwa index {input} 2> {log}"

rule make_genome_windows:
    input:
        f"Assembly/{bn}.fai"
    output: 
        f"Assembly/{bn}.bed"
    message: 
        "Creating BED intervals from {input}"
    shell: 
        "makewindows.py -i {input} -w 10000 -o {output}"

rule align:
    input:
        forward_reads = get_fq1,
        reverse_reads = get_fq2,
        genome 		  = f"Assembly/{bn}",
        genome_idx 	  = multiext(f"Assembly/{bn}", ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
    output:  
        temp(outdir + "/{sample}/{sample}.sort.bam"),
    log:
        outdir + "/logs/{sample}.log"
    message:
        "Aligning sequences: {wildcards.sample}"
    benchmark:
        "Benchmark/Mapping/bwa/align.{sample}.txt"
    params: 
        quality = config["quality"],
        tmpdir = lambda wc: outdir + "/" + d[wc.sample],
        extra   = extra
    threads:
        8
    shell:
        """
        #Align/bwa/{wildcards.sample}
        BWA_THREADS=$(( {threads} - 2 ))
        bwa mem -C -t $BWA_THREADS {params.extra} -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" {input.genome} {input.forward_reads} {input.reverse_reads} 2> {log} |
        samtools view -h -q {params.quality} | 
        samtools sort -T {params.tmpdir} --reference {input.genome} -O bam -l 0 -m 4G -o {output} 2> /dev/null
        rm -rf {params.tmpdir}
        """

rule mark_duplicates:
    input:
        lambda wc: outdir + "/" + d[wc.sample] + ".sort.bam"
    output:
        bam = temp(outdir + "/{sample}/{sample}.markdup.bam"),
        bai = temp(outdir + "/{sample}/{sample}.markdup.bam.bai")
    log:
        outdir + "/logs/{sample}.markdup.log"
    message:
        f"Marking duplicates: " + "{wildcards.sample}"
    benchmark:
        "Benchmark/Mapping/bwa/markdup.{sample}.txt"
    threads: 
        4
    shell:
        "sambamba markdup -t {threads} -l 0 {input} {output.bam} 2> {log}"

rule clip_overlap:
    input:
        bam = outdir + "/{sample}/{sample}.markdup.bam",
        bai = outdir + "/{sample}/{sample}.markdup.bam.bai"
    output:
        bam = outdir + "/{sample}.bam",
        bai = outdir + "/{sample}.bam.bai"
    log:
        outdir + "/logs/{sample}.clipOverlap.log"
    message:
        "Clipping alignment overlaps: {wildcards.sample}"
    shell:
        """
        bam clipOverlap --in {input.bam} --out {output.bam} --stats --noPhoneHome > {log} 2>&1
        samtools index {output.bam}
        """

rule alignment_coverage:
    input: 
        bed = f"Assembly/{bn}.bed",
        bam = outdir + "/{sample}.bam"
    output: 
        outdir + "/stats/coverage/data/{sample}.cov.gz"
    message:
        "Calculating genomic coverage: {wildcards.sample}"
    params:
        lambda wc: d[wc.sample]
    threads: 
        2
    shell:
        "samtools bedcov -c {input} | gzip > {output}"

rule coverage_report:
    input:
        outdir + "/stats/coverage/data/{sample}.cov.gz"
    output:
        outdir + "/stats/coverage/{sample}.cov.html"
    message:
        "Summarizing alignment coverage: {wildcards.sample}"
    script:
        "reportBwaGencov.Rmd"

rule alignment_bxstats:
    input:
        bam = outdir + "/{sample}.bam",
        bai = outdir + "/{sample}.bam.bai"
    output: 
        outdir + "/stats/BXstats/data/{sample}.bxstats.gz"
    message:
        "Calculating barcode alignment statistics: {wildcards.sample}"
    params:
        sample = lambda wc: d[wc.sample],
    shell:
        "bxStats.py {input.bam} | gzip > {output}"

rule bx_stats_report:
    input:
        outdir + "/stats/BXstats/data/{sample}.bxstats.gz"
    output:	
        outdir + "/stats/BXstats/{sample}.bxstats.html"
    message: 
        "Generating summary of barcode alignment: {wildcards.sample}"
    params:
        sample = lambda wc: d[wc.sample]
    script:
        "reportBxStats.Rmd"
    
rule general_alignment_stats:
    input:
        bam      = outdir + "/{sample}.bam",
        bai      = outdir + "/{sample}.bam.bai"
    output: 
        stats    = temp(outdir + "/stats/samtools_stats/{sample}.stats"),
        flagstat = temp(outdir + "/stats/samtools_flagstat/{sample}.flagstat")
    message:
        "Calculating alignment stats: {wildcards.sample}"
    params:
        sample = lambda wc: d[wc.sample]
    benchmark:
        "Benchmark/Mapping/bwa/stats.{sample}.txt"
    shell:
        """
        samtools stats {input.bam} > {output.stats}
        samtools flagstat {input.bam} > {output.flagstat}
        """

rule samtools_reports:
    input: 
        expand(outdir + "/stats/samtools_{ext}/{sample}.{ext}", sample = samplenames, ext = ["stats", "flagstat"])
    output: 
        outdir + "/stats/bwa.stats.html",
    message:
        "Summarizing samtools stats and flagstat"
    shell:
        "multiqc Align/bwa/stats/samtools_stats Align/bwa/stats/samtools_flagstat --force --quiet --no-data-dir --filename {output} 2> /dev/null"

rule log_runtime:
    output:
        outdir + "/logs/harpy.align.log"
    message:
        "Creating record of relevant runtime parameters: {output}"
    params:
        quality = config["quality"],
        extra   = extra
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy align module ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write(f"The directory with sequences: {seq_dir}\n\n")
            _ = f.write("Sequencing were aligned with BWA using:\n")
            _ = f.write("\tbwa mem -C " + " ".join([str(i) for i in params]) + " -R \"@RG\\tID:SAMPLE\\tSM:SAMPLE\" genome forward_reads reverse_reads |\n")
            _ = f.write("\tsamtools view -h -q " + str(config["quality"]) + " |\n")
            _ = f.write("\tsamtools sort -T SAMPLE --reference genome -m 4G\n")
            _ = f.write("Duplicates in the alignments were marked using sambamba:\n")
            _ = f.write("\tsambamba markdup -l 0\n")