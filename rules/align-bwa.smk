import os
import re
import glob

outdir      = "Align/bwa"
seq_dir		= config["seq_directory"]
genomefile 	= config["genomefile"]
samplenames = config["samplenames"]
extra 		= config.get("extra", "") 
bn 			= os.path.basename(genomefile)

d = dict(zip(samplenames, samplenames))

def get_fq1(wildcards):
    # returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    # the list is just the single reverse file
    lst = glob.glob(seq_dir + "/" + wildcards.sample + "*")
    return lst[0]

def get_fq2(wildcards):
    # returns a list of fastq files for read 2 based on *wildcards.sample*, e.g.
    # the list is just the single reverse file
    lst = sorted(glob.glob(seq_dir + "/" + wildcards.sample + "*"))
    return lst[1]

rule all:
    input: 
        bam = expand(outdir + "/align/{sample}.bam", sample = samplenames),
        bai = expand(outdir + "/align/{sample}.bam.bai", sample = samplenames),
        covreport = expand(outdir + "/stats/coverage/{sample}.cov.html", sample = samplenames),
        bxreport = expand(outdir + "/stats/BXstats/{sample}.bxstats.html", sample = samplenames),
        statsreport = outdir + "/stats/bwa.stats.html",
        runlog = outdir + "/logs/harpy.align.log"
    message:
        "Finished aligning! Moving alignment files into the base Align/bwa directory."
    default_target: True
    run:
        for i,j in zip(input.bam, input.bai):
            if not os.path.islink(i):
                # yank out just the filename
                fname = os.path.basename(i)
                # move file into base path
                os.rename(i, f"{outdir}/{fname}")
                # preserve "original" in align folder as symlink
                target = Path(f"{outdir}/{fname}").absolute()
                _ = Path(i).symlink_to(target)
            if not os.path.islink(j):
                # same for .bai file
                fnamebai = os.path.basename(j)
                os.rename(j, f"{outdir}/{fnamebai}")
                targetbai = Path(f"{outdir}/{fnamebai}").absolute()
                _ = Path(j).symlink_to(targetbai)

rule genome_link:
    input:
        genomefile
    output: 
        f"Genome/{bn}"
    message: 
        "Symlinking {input}"
    shell: 
        "ln -sr {input} {output}"

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

rule genome_bwa_index:
    input: 
        f"Genome/{bn}"
    output: 
        multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    message:
        "Indexing {input}"
    log:
        f"Genome/{bn}.idx.log"
    shell: 
        "bwa index {input} 2> {log}"

rule genome_make_windows:
    input:
        f"Genome/{bn}.fai"
    output: 
        f"Genome/{bn}.bed"
    message: 
        "Creating BED intervals from {input}"
    shell: 
        "makewindows.py -i {input} -w 10000 -o {output}"

rule align:
    input:
        forward_reads = get_fq1,
        reverse_reads = get_fq2,
        genome 		  = f"Genome/{bn}",
        genome_idx 	  = multiext(f"Genome/{bn}", ".ann", ".bwt", ".fai", ".pac", ".sa", ".amb")
    output:  
        bam = temp(outdir + "/{sample}/{sample}.sort.bam"),
        bai = temp(outdir + "/{sample}/{sample}.sort.bam.bai")
    log:
        bwa     = outdir + "/logs/{sample}.bwa.align.log",
        bwasort = outdir + "/logs/{sample}.bwa.sort.log"
    message:
        "Aligning sequences: {wildcards.sample}"
    benchmark:
        "Benchmark/Mapping/bwa/align.{sample}.txt"
    params: 
        quality = config["quality"],
        tmpdir = lambda wc: outdir + "/." + d[wc.sample],
        extra   = extra
    threads:
        10
    shell:
        """
        BWA_THREADS=$(( {threads} - 2 ))
        bwa mem -C -t $BWA_THREADS {params.extra} -R \"@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\" {input.genome} {input.forward_reads} {input.reverse_reads} 2> {log.bwa} |
        samtools view -h -q {params.quality} | 
        samtools sort -T {params.tmpdir} --reference {input.genome} -O bam -l 0 -m 4G --write-index -o {output.bam}##idx##{output.bai} 2> {log.bwasort}
        rm -rf {params.tmpdir}
        """

rule alignment_bxstats:
    input:
        bam = outdir + "/{sample}/{sample}.sort.bam",
        bai = outdir + "/{sample}/{sample}.sort.bam.bai"
    output: 
        outdir + "/stats/BXstats/data/{sample}.bxstats.gz"
    message:
        "Calculating barcode alignment statistics: {wildcards.sample}"
    params:
        sample = lambda wc: d[wc.sample],
    shell:
        "bxStats.py {input.bam} | gzip > {output}"

rule bxstats_report:
    input:
        outdir + "/stats/BXstats/data/{sample}.bxstats.gz"
    output:	
        outdir + "/stats/BXstats/{sample}.bxstats.html"
    message: 
        "Generating summary of barcode alignment: {wildcards.sample}"
    script:
        "reportBxStats.Rmd"

rule mark_duplicates:
    input:
        lambda wc: outdir + "/{sample}/{sample}.sort.bam"
    output:
        bam = temp(outdir + "/{sample}/{sample}.markdup.bam"),
        bai = temp(outdir + "/{sample}/{sample}.markdup.bam.bai")
    log:
        outdir + "/logs/makrduplicates/{sample}.markdup.log"
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
        bam = outdir + "/align/{sample}.bam",
        bai = outdir + "/align/{sample}.bam.bai"
    log:
        outdir + "/logs/clipOverlap/{sample}.clipOverlap.log"
    message:
        "Clipping alignment overlaps: {wildcards.sample}"
    shell:
        """
        bam clipOverlap --in {input.bam} --out {output.bam} --stats --noPhoneHome > {log} 2>&1
        samtools index {output.bam}
        """

rule alignment_coverage:
    input: 
        bed = f"Genome/{bn}.bed",
        bam = outdir + "/align/{sample}.bam"
    output: 
        outdir + "/stats/coverage/data/{sample}.cov.gz"
    message:
        "Calculating genomic coverage: {wildcards.sample}"
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
    
rule general_alignment_stats:
    input:
        bam      = outdir + "/align/{sample}.bam",
        bai      = outdir + "/align/{sample}.bam.bai"
    output: 
        stats    = temp(outdir + "/stats/samtools_stats/{sample}.stats"),
        flagstat = temp(outdir + "/stats/samtools_flagstat/{sample}.flagstat")
    message:
        "Calculating alignment stats: {wildcards.sample}"
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
            _ = f.write("Overlaps were clipped using:\n")
            _ = f.write("\tbam clipOverlap --in file.bam --out outfile.bam --stats --noPhoneHome\n")
