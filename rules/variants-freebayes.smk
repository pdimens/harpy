import os
import sys
import gzip

bam_dir 	= config["seq_directory"]
genomefile 	= config["genomefile"]
groupings 	= config.get("groupings", None)
bn          = os.path.basename(genomefile)
ploidy 		= config["ploidy"]
samplenames = config["samplenames"]
extra 	    = config.get("extra", "") 
chunksize   = config["windowsize"]
intervals   = config["intervals"]
outdir      = "Variants/freebayes"
regions     = dict(zip(intervals, intervals))


if groupings:
    rule copy_groupings:
        input:
            groupings
        output:
            outdir + "/logs/sample.groups"
        message:
            "Logging {input}"
        run:
            with open(input[0], "r") as infile, open(output[0], "w") as outfile:
                _ = [outfile.write(i) for i in infile.readlines() if not i.lstrip().startswith("#")]

rule index_alignments:
    input:
        bam_dir + "/{sample}.bam"
    output:
        bam_dir + "/{sample}.bam.bai"
    message:
        "Indexing alignments: {wildcards.sample}"
    benchmark:
        ".Benchmark/Variants/mpileup/indexbam.{sample}.txt"
    shell:
        "sambamba index {input} {output} 2> /dev/null"

rule samplenames:
    output:
        outdir + "/logs/samples.names"
    message:
        "Creating list of sample names"
    run:
        with open(output[0], "w") as fout:
            for samplename in samplenames:
                _ = fout.write(samplename + "\n")	

rule bam_list:
    input: 
        bam = expand(bam_dir + "/{sample}.bam", sample = samplenames),
        bai = expand(bam_dir + "/{sample}.bam.bai", sample = samplenames)
    output:
        outdir + "/logs/samples.files"
    message:
        "Creating list of alignment files"
    run:
        with open(output[0], "w") as fout:
            for bamfile in input.bam:
                _ = fout.write(bamfile + "\n")

if groupings:
    rule call_variants_pop:
        input:
            bam = expand(bam_dir + "/{sample}.bam", sample = samplenames),
            bai = expand(bam_dir + "/{sample}.bam.bai", sample = samplenames),
            groupings = outdir + "/logs/sample.groups",
            ref     = f"Genome/{bn}",
            ref_idx = f"Genome/{bn}.fai",
            samples = outdir + "/logs/samples.files"
        output:
            bcf = temp(outdir + "/regions/{part}.bcf"),
            idx = temp(outdir + "/regions/{part}.bcf.csi")
        message:
            "Calling variants: {wildcards.part}"
        threads:
            2
        params:
            region = lambda wc: "-r " + regions[wc.part],
            ploidy = f"-p {ploidy}",
            extra = extra
        shell:
            """
            freebayes -f {input.ref} -L {input.samples} --populations {input.groupings} {params} | bcftools sort - -Ob --output {output.bcf} 2> /dev/null
            bcftools index {output.bcf}
            """
else:
    rule call_variants:
        input:
            bam = expand(bam_dir + "/{sample}.bam", sample = samplenames),
            bai = expand(bam_dir + "/{sample}.bam.bai", sample = samplenames),
            ref     = f"Genome/{bn}",
            ref_idx = f"Genome/{bn}.fai",
            samples = outdir + "/logs/samples.files"
        output:
            bcf = temp(outdir + "/regions/{part}.bcf"),
            idx = temp(outdir + "/regions/{part}.bcf.csi")
        message:
            "Calling variants: {wildcards.part}"
        threads:
            2
        params:
            region = lambda wc: "-r " + regions[wc.part],
            ploidy = f"-p {ploidy}",
            extra = extra
        shell:
            """
            freebayes -f {input.ref} -L {input.samples} {params} | bcftools sort - -Ob --output {output.bcf} 2> /dev/null
            bcftools index {output.bcf}
            """

rule concat_list:
    input:
        bcfs = expand(outdir + "/regions/{part}.bcf", part = intervals),
    output:
        outdir + "/logs/bcf.files"
    message:
        "Creating list of region-specific vcf files"
    run:
        with open(output[0], "w") as fout:
            for bcf in input.bcfs:
                _ = fout.write(f"{bcf}\n")   

rule merge_vcfs:
    input:
        bcfs = expand(outdir + "/regions/{part}.{ext}", part = intervals, ext = ["bcf", "bcf.csi"]),
        filelist = outdir + "/logs/bcf.files"
    output:
        bcf = outdir + "/variants.raw.bcf",
        idx = outdir + "/variants.raw.bcf.csi"
    message:
        "Combining vcfs into a single file"
    log:
        outdir + "/logs/concat.log"
    threads:
        50
    shell:  
        """
        #bcftools concat -f {input.filelist} --threads {threads} --naive -Ob --write-index > {output.bcf} 2> {log}
        bcftools concat -f {input.filelist} --threads {threads} --naive -Ob > {output.bcf} 2> {log}
        bcftools index --threads {threads} {output.bcf}
        """

rule normalize_bcf:
    input: 
        genome  = f"Genome/{bn}",
        ref_idx = f"Genome/{bn}.fai",
        bcf     = outdir + "/variants.raw.bcf"
    output:
        bcf     = outdir + "/variants.normalized.bcf",
        idx     = outdir + "/variants.normalized.bcf.csi",
    message: 
        "Normalizing the called variants"
    threads: 
        2
    shell:
        """
        #bcftools norm -d exact -f {input.genome} {input.bcf} | bcftools norm -m -any -N -Ob --write-index > {output.bcf}
        bcftools norm -d exact -f {input.genome} {input.bcf} | bcftools norm -m -any -N -Ob > {output.bcf}
        bcftools index --threads {threads} {output.bcf}        
        """

rule variants_stats:
    input:
        genome  = f"Genome/{bn}",
        ref_idx = f"Genome/{bn}.fai",
        bcf     = outdir + "/variants.{type}.bcf",
        idx     = outdir + "/variants.{type}.bcf.csi"
    output:
        outdir + "/stats/variants.{type}.stats",
    message:
        "Calculating variant stats: variants.{wildcards.type}.bcf"
    shell:
        """bcftools stats -s "-" --fasta-ref {input.genome} {input.bcf} > {output}"""

rule bcfreport:
    input:
        outdir + "/stats/variants.{type}.stats"
    output:
        outdir + "/stats/variants.{type}.html"
    message:
        "Generating bcftools report: variants.{wildcards.type}.bcf"
    script:
        "reportBcftools.Rmd"

rule log_runtime:
    output:
        outdir + "/logs/harpy.variants.log"
    message:
        "Creating record of relevant runtime parameters: {output}"
    params:
        ploidy = f"-p {ploidy}",
        populations = '' if groupings is None else f"--populations {groupings}",
        extra = extra
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy variants snp module ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write(f"The directory with alignments: {bam_dir}\n")
            _ = f.write(f"Size of intervals to split genome for variant calling: {chunksize}\n")
            _ = f.write("The freebayes parameters:\n")
            _ = f.write("    freebayes -f GENOME -L samples.list -r REGION " + " ".join(params) + " | bcftools sort -\n")
            _ = f.write("The variants identified in the intervals were merged into the final variant file using:\n")
            _ = f.write("    bcftools concat -f vcf.list -a --remove-duplicates\n")
            _ = f.write("The variants were normalized using:\n")
            _ = f.write("    bcftools norm -d exact | bcftools norm -m -any -N -Ob\n")

rule all:
    default_target: True
    input: 
        outdir + "/logs/harpy.variants.log",
        expand(outdir + "/variants.{file}.bcf",        file = ["raw", "normalized"]),
        expand(outdir + "/stats/variants.{file}.html", file = ["raw", "normalized"])
    message:
        "Variant calling is complete!"