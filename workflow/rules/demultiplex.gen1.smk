import os
import re
import sys
from rich.panel import Panel
from rich import print as rprint

infile = config["infile"]
samplefile = config["samplefile"]
skipreports = config["skipreports"]
bn = os.path.basename(infile)
fq_extension = re.search(r"(?:\_00[0-9])*\.f(.*?)q(?:\.gz)?$", infile, re.IGNORECASE).group(0)
inprefix = config["infile_prefix"]
inprefixfull = re.sub(r"[\_\.][IR][12]?(?:\_00[0-9])*\.f(?:ast)?q(?:\.gz)?$", "", infile)
infiles = [f"{inprefixfull}_{i}{fq_extension}" for i in ["I1", "I2","R1","R2"]]
indir = os.path.dirname(infile)
outdir = f"Demultiplex/{inprefix}/"

def barcodedict(smpl):
    d = dict()
    with open(smpl, "r") as f:
        for i in f.readlines():
            # a casual way to ignore empty lines or lines with >2 fields
            try:
                smpl, bc = i.split()
                d[smpl] = bc
            except:
                continue
    return d

samples = barcodedict(samplefile)
samplenames = [i for i in samples.keys()]

onerror:
    print("")
    rprint(
        Panel(
            f"The workflow has terminated due to an error. See the log file below for more details.",
            title = "[bold]harpy demultiplex",
            title_align = "left",
            border_style = "red"
            ),
        file = sys.stderr
    )

onsuccess:
    print("")
    rprint(
        Panel(
            f"The workflow has finished successfully! Find the results in [bold]{outdir}[/bold]",
            title = "[bold]harpy demultiplex",
            title_align = "left",
            border_style = "green"
            ),
        file = sys.stderr
    )

rule link_files:
    input:
        indir + "/" + inprefix + "_{part}" + fq_extension
    output:
        temp(outdir + "DATA_{part}_001.fastq.gz")
    message:
        f"Linking {inprefix}" + "{wildcards.part} to output directory"
    shell:
        "ln -sr {input} {output}"

rule bx_files:
    output:
        temp(expand(outdir + "BC_{letter}.txt", letter = ["A","C","B","D"]))
    params:
        outdr = outdir
    message:
        "Creating the Gen I barcode files necessary for barcode demultiplexing"
    shell:
        """
        cd {params}
        BC_files.py
        """

rule demux_bx:
    input:
        expand(outdir + "DATA_{IR}{ext}_001.fastq.gz", IR = ["R","I"], ext = [1,2]),
        expand(outdir + "BC_{letter}.txt", letter = ["A","C","B","D"])
    output:
        temp(expand(outdir + inprefix + "_R{ext}_001.fastq.gz", ext = [1,2]))
    params:
        outdr = outdir,
        outprfx = inprefix,
        logdir = outdir +"logs/.QC"
    message:
        "Moving barcodes into read headers"
    shell:
        """
        mkdir -p {params.logdir}
        cd {params.outdr}
        demuxGen1 DATA_ {params.outprfx}
        mv {params.outprfx}*BC.log logs
        """

rule split_samples_fw:
    input:
        f"{outdir}{inprefix}_R1_001.fastq.gz"
    output:
        outdir + "{sample}.F.fq.gz"
    params:
        c_barcode = lambda wc: samples[wc.get("sample")]
    message:
        "Extracting forward reads:\n sample: {wildcards.sample}\n barcode: {params}"
    shell:
        """
        ( zgrep -A3 "A..{params}B..D" {input} | grep -v "^--$" | gzip -q > {output} ) || touch {output}
        """

rule split_samples_rv:
    input:
        f"{outdir}{inprefix}_R2_001.fastq.gz"
    output:
        outdir + "{sample}.R.fq.gz"
    params:
        c_barcode = lambda wc: samples[wc.get("sample")]
    message:
        "Extracting reverse reads:\n sample: {wildcards.sample}\n barcode: {params}"
    shell:
        """
        ( zgrep -A3 "A..{params}B..D" {input} | grep -v "^--$" | gzip -q > {output} ) || touch {output}
        """

rule fastqc_F:
    input:
        outdir + "{sample}.F.fq.gz"
    output: 
        temp(outdir + "logs/.QC/{sample}_F/fastqc_data.txt")
    params:
        lambda wc: outdir + "logs/.QC/" + wc.get("sample") + "_F"
    threads:
        1
    conda:
        os.getcwd() + "/.harpy_envs/qc.yaml"
    message:
        "Performing quality assessment: {wildcards.sample}.F.fq.gz"
    shell:
        """
        mkdir -p {params}
        if [ -z $(gzip -cd {input} | head -c1) ]; then
            echo "##Falco	1.2.1" > {output}
            echo ">>Basic Statistics	fail" >> {output}
            echo "#Measure	Value" >> {output}
            echo "Filename	{wildcards.sample}.F.fq.gz" >> {output}
            echo "File type	Conventional base calls" >> {output}
            echo "Encoding	Sanger / Illumina 1.9" >> {output}
            echo "Total Sequences	0" >> {output}
            echo "Sequences flagged as poor quality	0" >> {output}
            echo "Sequence length	0" >> {output}
            echo "%GC	0" >> {output}
            echo ">>END_MODULE" >> {output}
        else
            falco -q --threads {threads} -skip-report -skip-summary -o {params} {input}
        fi
        """

rule fastqc_R:
    input:
        outdir + "{sample}.R.fq.gz"
    output: 
        temp(outdir + "logs/.QC/{sample}_R/fastqc_data.txt")
    params:
        lambda wc: outdir + "logs/.QC/" + wc.get("sample") + "_R"
    threads:
        1
    conda:
        os.getcwd() + "/.harpy_envs/qc.yaml"
    message:
        "Performing quality assessment: {wildcards.sample}.R.fq.gz"
    shell:
        """
        mkdir -p {params}
        if [ -z $(gzip -cd {input} | head -c1) ]; then
            echo "##Falco	1.2.1" > {output}
            echo ">>Basic Statistics	fail" >> {output}
            echo "#Measure	Value" >> {output}
            echo "Filename	{wildcards.sample}.F.fq.gz" >> {output}
            echo "File type	Conventional base calls" >> {output}
            echo "Encoding	Sanger / Illumina 1.9" >> {output}
            echo "Total Sequences	0" >> {output}
            echo "Sequences flagged as poor quality	0" >> {output}
            echo "Sequence length	0" >> {output}
            echo "%GC	0" >> {output}
            echo ">>END_MODULE" >> {output}
        else
            falco -q --threads {threads} -skip-report -skip-summary -o {params} {input}
        fi
        """

rule qc_report:
    input:
        expand(outdir + "logs/.QC/{sample}_{FR}/fastqc_data.txt", sample = samplenames, FR = ["F","R"])
    output:
        outdir + "reports/demultiplex.QC.html"
    params:
        outdir + "logs/.QC"
    conda:
        os.getcwd() + "/.harpy_envs/qc.yaml"
    message:
        "Creating final demultiplexing QC report"
    shell:
        """
        multiqc {params} --force --quiet --title "QC on Demultiplexed Samples" --comment "This report aggregates the QC results created by falco." --no-data-dir --filename {output} 2> /dev/null
        """

rule log_runtime:
    output:
        outdir + "workflow/demultiplex.workflow.summary"
    message:
        "Creating record of relevant runtime parameters: {output}"
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy demultiplex module ran using these parameters:\n\n")
            _ = f.write("Haplotag technology: Generation I\n")
            _ = f.write(f"The multiplexed input file: {infile}\n")
            _ = f.write(f"The associated files inferred from {infile}:\n")
            _ = f.write("    " + "\n    ".join(infiles) + "\n")
            _ = f.write("Barcodes were moved into the read headers using the command:\n")
            _ = f.write(f"    demuxGen1 DATA_ {inprefix}\n")
            _ = f.write(f"The delimited file associating CXX barcodes with samplenames: {samplefile}\n")
            _ = f.write(f"QC checks were performed on demultiplexed FASTQ files using:\n")
            _ = f.write(f"    falco -skip-report -skip-summary input.fq.gz\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")

# conditionally add the reports to the output
results = list()
results.append(expand(outdir + "{sample}.F.fq.gz", sample = samplenames))
results.append(expand(outdir + "{sample}.R.fq.gz", sample = samplenames))
results.append(outdir + "workflow/demultiplex.workflow.summary")

if not skipreports:
    results.append(outdir + "reports/demultiplex.QC.html")

rule all:
    default_target: True
    input:
        results
    message:
        "Checking for expected workflow output"
