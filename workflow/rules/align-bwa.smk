import os
import re
import glob
import sys
from rich.panel import Panel
from rich import print as rprint

outdir      = config["output_directory"]
seq_dir		= config["seq_directory"]
genomefile 	= config["genomefile"]
samplenames = config["samplenames"]
molecule_distance = config["molecule_distance"]
extra 		= config.get("extra", "") 
bn 			= os.path.basename(genomefile)
genome_zip  = True if bn.lower().endswith(".gz") else False
bn_idx      = f"{bn}.gzi" if genome_zip else f"{bn}.fai"
skipreports = config["skipreports"]

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
        log:
            f"Genome/{bn}.faidx.gzi.log"
        message:
            "Indexing {input}"
        shell: 
            "samtools faidx --gzi-idx {output.gzi} --fai-idx {output.fai} {input} 2> {log}"
else:
    rule genome_faidx:
        input: 
            f"Genome/{bn}"
        output: 
            f"Genome/{bn}.fai"
        log:
            f"Genome/{bn}.faidx.log"
        message:
            "Indexing {input}"
        shell:
            "samtools faidx --fai-idx {output} {input} 2> {log}"

rule genome_bwa_index:
    input: 
        f"Genome/{bn}"
    output: 
        multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    log:
        f"Genome/{bn}.idx.log"
    conda:
        os.getcwd() + "/.harpy_envs/align.yaml"
    message:
        "Indexing {input}"
    shell: 
        "bwa index {input} 2> {log}"

rule genome_make_windows:
    input:
        f"Genome/{bn}"
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
        genome_idx 	  = multiext(f"Genome/{bn}", ".ann", ".bwt", ".pac", ".sa", ".amb")
    output:  
        pipe(outdir + "/samples/{sample}/{sample}.raw.sam")
    log:
        outdir + "/logs/{sample}.bwa.log"
    params: 
        samps = lambda wc: d[wc.get("sample")],
        extra = extra
    benchmark:
        ".Benchmark/Mapping/bwa/align.{sample}.txt"
    threads:
        min(10, workflow.cores) - 2
    conda:
        os.getcwd() + "/.harpy_envs/align.yaml"
    message:
        "Aligning sequences: {wildcards.sample}"
    benchmark:
        ".Benchmark/Mapping/bwa/align.{sample}.txt"
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
        samtools view -h -F 4 -q {params.quality} | 
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
        mdist = molecule_distance
    shell:
        "bxStats.py -c {params.mdist} {input.bam} | gzip > {output}"

rule bxstats_report:
    input:
        outdir + "/reports/BXstats/data/{sample}.bxstats.gz"
    output:	
        outdir + "/reports/BXstats/{sample}.bxstats.html"
    params:
        molecule_distance
    conda:
        os.getcwd() + "/.harpy_envs/r-env.yaml"
    message: 
        "Generating summary of barcode alignment: {wildcards.sample}"
    params:
        molecule_distance
    script:
        "report/BxStats.Rmd"

rule mark_duplicates:
    input:
        outdir + "/samples/{sample}/{sample}.sort.bam"
    output:
        bam = temp(outdir + "/{sample}/markdup/{sample}.markdup.bam"),
        bai = temp(outdir + "/{sample}/markdup/{sample}.markdup.bam.bai")
    log:
        outdir + "/logs/makrduplicates/{sample}.markdup.log"
    message:
        f"Marking duplicates: " + "{wildcards.sample}"
    benchmark:
        ".Benchmark/Mapping/bwa/markdup.{sample}.txt"
    threads: 
        4
    conda:
        os.getcwd() + "/harpyenvs/align.yaml"
    message:
        "Marking duplicates: " + "{wildcards.sample}"
    shell:
        "sambamba markdup -t {threads} -l 0 {input} {output.bam} 2> {log}"

rule assign_molecules:
    input:
        bam = outdir + "/{sample}/markdup/{sample}.markdup.bam",
        bai = outdir + "/{sample}/markdup/{sample}.markdup.bam.bai"
    output:
        bam = outdir + "/align/{sample}.bam",
        bai = outdir + "/align/{sample}.bam.bai"
    message:
        "Assigning barcodes to molecules: {wildcards.sample}"
    params:
        molecule_distance
    shell:
        "assignMI.py -c {params} -i {input.bam} -o {output.bam}"

rule alignment_coverage:
    input: 
        bed = f"Genome/{bn}.bed",
        bam = outdir + "/{sample}.bam"
    output: 
        outdir + "/reports/coverage/data/{sample}.cov.gz"
    threads: 
        2
    message:
        "Calculating genomic coverage: {wildcards.sample}"
    shell:
        "samtools bedcov -c {input} | gzip > {output}"

rule coverage_report:
    input:
        outdir + "/reports/coverage/data/{sample}.cov.gz"
    output:
        outdir + "/reports/coverage/{sample}.cov.html"
    conda:
        os.getcwd() + "/.harpy_envs/r-env.yaml"
    message:
        "Summarizing alignment coverage: {wildcards.sample}"
    script:
        "report/Gencov.Rmd"
    
rule general_alignment_stats:
    input:
        bam      = outdir + "/{sample}.bam",
        bai      = outdir + "/{sample}.bam.bai"
    output: 
        stats    = temp(outdir + "/reports/samtools_stats/{sample}.stats"),
        flagstat = temp(outdir + "/reports/samtools_flagstat/{sample}.flagstat")
    message:
        "Calculating alignment stats: {wildcards.sample}"
    benchmark:
        ".Benchmark/Mapping/bwa/stats.{sample}.txt"
    shell:
        """
        samtools stats {input.bam} > {output.stats}
        samtools flagstat {input.bam} > {output.flagstat}
        """

rule samtools_reports:
    input: 
        expand(outdir + "/reports/samtools_{ext}/{sample}.{ext}", sample = samplenames, ext = ["stats", "flagstat"])
    output: 
        outdir + "/reports/bwa.stats.html"
    params:
        outdir
    conda:
        os.getcwd() + "/.harpy_envs/qc.yaml"
    message:
        "Summarizing samtools stats and flagstat"
    shell:
        """
        multiqc Align/bwa/stats/samtools_stats Align/bwa/stats/samtools_flagstat --force --quiet --title "Basic Alignment Statistics" --comment "This report aggregates samtools stats and samtools flagstats results for all alignments." --no-data-dir --filename {output} 2> /dev/null
        """

rule log_runtime:
    output:
        outdir + "/workflow/align.workflow.summary"
    params:
        quality = config["quality"],
        extra   = extra
    message:
        "Creating record of relevant runtime parameters: {output}"
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy align module ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write(f"The directory with sequences: {seq_dir}\n\n")
            _ = f.write("Sequencing were aligned with BWA using:\n")
            _ = f.write("    bwa mem -C " + " ".join([str(i) for i in params]) + " -R \"@RG\\tID:SAMPLE\\tSM:SAMPLE\" genome forward_reads reverse_reads |\n")
            _ = f.write("    samtools view -h -F 4 -q " + str(config["quality"]) + " |\n")
            _ = f.write("    samtools sort -T SAMPLE --reference genome -m 4G\n")
            _ = f.write("Duplicates in the alignments were marked using sambamba:\n")
            _ = f.write("    sambamba markdup -l 0\n")
            _ = f.write("Overlaps were clipped using:\n")
            _ = f.write("    bam clipOverlap --in file.bam --out outfile.bam --stats --noPhoneHome\n")

rule movelinks:
    default_target: True
    input: 
        bam = expand(outdir + "/align/{sample}.bam", sample = samplenames),
        bai = expand(outdir + "/align/{sample}.bam.bai", sample = samplenames),
        covreport = expand(outdir + "/stats/coverage/{sample}.cov.html", sample = samplenames),
        bxreport = expand(outdir + "/stats/BXstats/{sample}.bxstats.html", sample = samplenames),
        statsreport = outdir + "/stats/bwa.stats.html",
        runlog = outdir + "/logs/harpy.align.log"
    message:
        "Finished aligning! Moving alignment files into the base Align/bwa directory."
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