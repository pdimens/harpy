containerized: "docker://pdimens/harpy:latest"

from rich import print as rprint
from rich.panel import Panel
import os
import sys
import gzip

bam_dir 	= config["seq_directory"]
envdir      = os.getcwd() + "/.harpy_envs"
genomefile 	= config["genomefile"]
groupings 	= config.get("groupings", [])
bn          = os.path.basename(genomefile)
if bn.lower().endswith(".gz"):
    genome_zip  = True
    bn = bn[:-3]
else:
    genome_zip  = False

ploidy 		= config["ploidy"]
samplenames = config["samplenames"]
extra 	    = config.get("extra", "") 
regioninput = config["regions"]
regiontype  = config["regiontype"]
windowsize  = config.get("windowsize", None)
outdir      = config["output_directory"]
skipreports = config["skipreports"]

if regiontype == "region":
    intervals = [regioninput]
    regions = {f"{regioninput}" : f"{regioninput}"}
else:
    with open(regioninput, "r") as reg_in:
        intervals = set()
        while True:
            line = reg_in.readline()
            if not line:
                break
            cont,startpos,endpos = line.split()
            intervals.add(f"{cont}:{startpos}-{endpos}")
    regions = dict(zip(intervals, intervals))

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

onerror:
    print("")
    rprint(
        Panel(
            f"The workflow has terminated due to an error. See the log file below for more details.",
            title = "[bold]harpy snp freebayes",
            title_align = "left",
            border_style = "red"
            ),
        file = sys.stderr
    )

onsuccess:
    print("")
    rprint(
        Panel(
            f"The workflow has finished successfully! Find the results in [bold]{outdir}/[/bold]",
            title = "[bold]harpy snp freebayes",
            title_align = "left",
            border_style = "green"
            ),
        file = sys.stderr
    )

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

rule genome_link:
    input:
        genomefile
    output: 
        f"Genome/{bn}"
    container:
        None
    message: 
        "Symlinking {input}"
    shell: 
        """
        if (file {input} | grep -q compressed ) ;then
            # is regular gzipped, needs decompression
            gzip -d -c {input} > {output}
        elif (file {input} | grep -q BGZF ); then
            # is bgzipped, needs decompression
            gzip -d -c {input} > {output}
        else
            # isn't compressed, just linked
            ln -sr {input} {output}
        fi
        """

rule genome_faidx:
    input: 
        f"Genome/{bn}"
    output: 
        f"Genome/{bn}.fai"
    log:
        f"Genome/{bn}.faidx.log"
    container:
        None
    message:
        "Indexing {input}"
    shell:
        "samtools faidx --fai-idx {output} {input} 2> {log}"

rule index_alignments:
    input:
        bam_dir + "/{sample}.bam"
    output:
        bam_dir + "/{sample}.bam.bai"
    container:
        None
    message:
        "Indexing alignments: {wildcards.sample}"
    shell:
        "samtools index {input} {output} 2> /dev/null"

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
        bam = collect(bam_dir + "/{sample}.bam", sample = samplenames),
        bai = collect(bam_dir + "/{sample}.bam.bai", sample = samplenames)
    output:
        outdir + "/logs/samples.files"
    message:
        "Creating list of alignment files"
    run:
        with open(output[0], "w") as fout:
            for bamfile in input.bam:
                _ = fout.write(bamfile + "\n")

rule call_variants:
    input:
        bam = collect(bam_dir + "/{sample}.bam", sample = samplenames),
        bai = collect(bam_dir + "/{sample}.bam.bai", sample = samplenames),
        groupfile = outdir + "/logs/sample.groups" if groupings else [],
        ref     = f"Genome/{bn}",
        ref_idx = f"Genome/{bn}.fai",
        samples = outdir + "/logs/samples.files"
    output:
        pipe(outdir + "/regions/{part}.vcf")
    params:
        region = lambda wc: "-r " + regions[wc.part],
        ploidy = f"-p {ploidy}",
        populations = "--populations " + rules.copy_groupings.output[0] if groupings else "",
        extra = extra
    conda:
        f"{envdir}/snp.yaml"
    message:
        "Calling variants: {wildcards.part}"
    shell:
        "freebayes -f {input.ref} -L {input.samples} {params} > {output}"

rule sort_variants:
    input:
        outdir + "/regions/{part}.vcf"
    output:
        bcf = temp(outdir + "/regions/{part}.bcf"),
        idx = temp(outdir + "/regions/{part}.bcf.csi")
    container:
        None
    message:
        "Sorting: {wildcards.part}"
    shell:
        "bcftools sort -Ob --write-index --output {output.bcf} {input} 2> /dev/null"

rule concat_list:
    input:
        bcfs = collect(outdir + "/regions/{part}.bcf", part = intervals),
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
        bcfs = collect(outdir + "/regions/{part}.{ext}", part = intervals, ext = ["bcf", "bcf.csi"]),
        filelist = outdir + "/logs/bcf.files"
    output:
        temp(outdir + "/variants.raw.unsort.bcf")
    log:
        outdir + "/logs/concat.log"
    threads:
        workflow.cores
    container:
        None
    message:
        "Combining vcfs into a single file"
    shell:  
        "bcftools concat -f {input.filelist} --threads {threads} --naive -Ob -o {output} 2> {log}"

rule sort_vcf:
    input:
        outdir + "/variants.raw.unsort.bcf"
    output:
        bcf = outdir + "/variants.raw.bcf",
        csi = outdir + "/variants.raw.bcf.csi"
    container:
        None
    message:
        "Sorting and indexing final variants"
    shell:
        "bcftools sort --write-index -Ob -o {output.bcf} {input} 2> /dev/null"


#rule normalize_bcf:
#    input: 
#        genome  = f"Genome/{bn}",
#        ref_idx = f"Genome/{bn}.fai",
#        bcf     = outdir + "/variants.raw.bcf"
#    output:
#        bcf     = outdir + "/variants.normalized.bcf",
#        idx     = outdir + "/variants.normalized.bcf.csi"
#    log:
#        outdir + "/logs/normalize.log"
#    threads: 
#        2
#    message: 
#        "Normalizing the called variants"
#    shell:
#        """
#        bcftools norm -d exact -f {input.genome} {input.bcf} 2> {log}.tmp1 | 
#            bcftools norm -m -any -N -Ob --write-index -o {output.bcf} 2> {log}.tmp2
#        cat {log}.tmp1 {log}.tmp2 > {log} && rm {log}.tmp1 {log}.tmp2    
#        """

rule variants_stats:
    input:
        genome  = f"Genome/{bn}",
        ref_idx = f"Genome/{bn}.fai",
        bcf     = outdir + "/variants.{type}.bcf",
        idx     = outdir + "/variants.{type}.bcf.csi"
    output:
        outdir + "/reports/variants.{type}.stats",
    container:
        None
    message:
        "Calculating variant stats: variants.{wildcards.type}.bcf"
    shell:
        """
        bcftools stats -s "-" --fasta-ref {input.genome} {input.bcf} > {output} 2> /dev/null
        """

rule bcf_report:
    input:
        outdir + "/reports/variants.{type}.stats"
    output:
        outdir + "/reports/variants.{type}.html"
    conda:
        f"{envdir}/r.yaml"
    message:
        "Generating bcftools report: variants.{wildcards.type}.bcf"
    script:
        "report/BcftoolsStats.Rmd"

rule log_workflow:
    default_target: True
    input:
        vcf = collect(outdir + "/variants.{file}.bcf", file = ["raw"]),
        reports = collect(outdir + "/reports/variants.{file}.html", file = ["raw"]) if not skipreports else []
    message:
        "Summarizing the workflow: {output}"
    params:
        ploidy = f"-p {ploidy}",
        populations = f"--populations {groupings}" if groupings else '',
        extra = extra
    run:
        with open(outdir + "/workflow/snp.freebayes.summary", "w") as f:
            _ = f.write("The harpy variants snp module ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write(f"The directory with alignments: {bam_dir}\n")
            _ = f.write(f"Size of intervals to split genome for variant calling: {windowsize}\n")
            _ = f.write("The freebayes parameters:\n")
            _ = f.write("    freebayes -f GENOME -L samples.list -r REGION " + " ".join(params) + " | bcftools sort -\n")
            _ = f.write("The variants identified in the intervals were merged into the final variant file using:\n")
            _ = f.write("    bcftools concat -f vcf.list -a --remove-duplicates\n")
            #_ = f.write("The variants were normalized using:\n")
            #_ = f.write("    bcftools norm -d exact | bcftools norm -m -any -N -Ob\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")