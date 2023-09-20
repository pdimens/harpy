import os
import re
import sys

infile = config["infile"]
samplefile = config["samplefile"]
#paramspace  = Paramspace(pd.read_csv(samplefile, sep="\t", header = None, names = ["barcode", "sample"]), param_sep = "", filename_params="*")
bn = os.path.basename(infile)
fq_extension = re.search(r"(?:\_00[0-9])*\.f(.*?)q(?:\.gz)?$", infile, re.IGNORECASE).group(0)
inprefix     = re.sub(r"[\_\.][IR][12]?(?:\_00[0-9])*\.f(?:ast)?q(?:\.gz)?$", "", bn)
inprefixfull = re.sub(r"[\_\.][IR][12]?(?:\_00[0-9])*\.f(?:ast)?q(?:\.gz)?$", "", infile)
infiles = [f"{inprefixfull}_{i}{fq_extension}" for i in ["I1", "I2","R1","R2"]]
indir = os.path.dirname(infile)
outdir = f"Demultiplex/{inprefix}/"


# mv functions to harpy executable?
## depends on what other demux methods look like
def checkfiles(prefix, prefixfull, ext):
    filelist = []
    printerr = False
    for i in ["I1", "I2","R1","R2"]:
        chkfile = f"{prefixfull}_{i}{ext}"
        TF = os.path.exists(chkfile)
        printerr = True if not TF else printerr
        symbol = " " if TF else "X"
        filelist.append(f"\033[91m{symbol}\033[0m  {prefix}_{i}{ext}")
    if printerr:
        print(f"\n\033[91mError\033[0m: Not all necessary files with prefix \033[1m{prefix}\033[0m present")
        _ = [print(i, file = sys.stderr) for i in filelist]
        exit(1)


checkfiles(inprefix, inprefixfull, fq_extension)

def get_samplenames(smpl):
    d = dict()
    with open(smpl, "r") as f:
        #rows = [i.split("\t")[0] for i in f.readlines()]
        for i in f.readlines():
            # a casual way to ignore empty lines or lines with >2 fields
            try:
                smpl, bc = i.split()
                d[smpl] = bc
            except:
                continue
    return d

samples = get_samplenames(samplefile)
samplenames = [i for i in samples.keys()]

rule link_files:
    input:
        indir + "/" + inprefix + "_{part}" + fq_extension
    output:
        temp(outdir + "DATA_{part}_001.fastq.gz")
    message:
        f"Linking {inprefix}" + "{wildcards.part} to output directory"
    shell:
        """
        ln -sr {input} {output}
        """

rule bx_files:
    output:
        temp(expand(outdir + "BC_{letter}.txt", letter = ["A","C","B","D"]))
    message:
        "Creating the Gen I barcode files necessary for barcode demultiplexing"
    params:
        outdr = outdir
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
    message:
        "Moving barcodes into read headers"
    params:
        outdr = outdir,
        outprfx = inprefix,
        logdir = outdir +"logs/.QC"
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
    message:
        "Demultiplexing forward reads: {wildcards.sample}"
    params:
        c_barcode = lambda wc: samples[wc.get("sample")]
    shell:
        """
        ( zgrep -A3 "A..{params}B..D" {input} | grep -v "^--$" | gzip -q > {output} ) || touch {output}
        """

rule split_samples_rv:
    input:
        f"{outdir}{inprefix}_R2_001.fastq.gz"
    output:
        outdir + "{sample}.R.fq.gz"
    message:
        "Demultiplexing reverse reads: {wildcards.sample}"
    params:
        c_barcode = lambda wc: samples[wc.get("sample")]
    shell:
        """
        ( zgrep -A3 "A..{params}B..D" {input} | grep -v "^--$" | gzip -q > {output} ) || touch {output}
        """

rule fastqc_F:
    input:
        outdir + "{sample}.F.fq.gz"
    output: 
        temp(outdir + "logs/.QC/{sample}_F/fastqc_data.txt")
    message:
        "Performing quality assessment: {wildcards.sample}.F.fq.gz"
    params:
        lambda wc: outdir + "logs/.QC/" + wc.get("sample") + "_F"
    threads:
        1
    shell:
        """
        mkdir -p {params}
        if [ -z $(gzip -cd {input} | head -c1) ]; then
            touch {output}
        else
            falco -q --threads {threads} -skip-report -skip-summary -o {params} {input}
        fi
        """

rule fastqc_R:
    input:
        outdir + "{sample}.R.fq.gz"
    output: 
        temp(outdir + "logs/.QC/{sample}_R/fastqc_data.txt")
    message:
        "Performing quality assessment: {wildcards.sample}.R.fq.gz"
    params:
        lambda wc: outdir + "logs/.QC/" + wc.get("sample") + "_R"
    threads:
        2
    shell:
        """
        mkdir -p {params}
        if [ -z $(gzip -cd {input} | head -c1) ]; then
            touch {output}
        else
            falco -q --threads {threads} -skip-report -skip-summary -o {params} {input}
        fi
        """

rule qc_report:
    input:
        expand(outdir + "logs/.QC/{sample}_{FR}/fastqc_data.txt", sample = samplenames, FR = ["F","R"])
    output:
        outdir + "logs/demultiplex.QC.html"
    message:
        "Creating final demultiplexing QC report"
    params:
        outdir + "logs/.QC"
    shell:
        "multiqc {params} --force --quiet --no-data-dir --filename {output} 2> /dev/null"

rule log_runtime:
    output:
        outdir + "logs/harpy.demultiplex.log"
    message:
        "Creating record of relevant runtime parameters: {output}"
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy demultiplex module ran using these parameters:\n\n")
            _ = f.write("Haplotag technology: Generation I\n")
            _ = f.write(f"The multiplexed input file: {infile}\n")
            _ = f.write(f"The inferred files associated with {infile}:\n")
            _ = f.write("    " + "\n    ".join(infiles) + "\n")
            _ = f.write("Barcodes were moved into the read headers using the command:\n")
            _ = f.write(f"    demuxGen1 DATA_ {inprefix}\n")
            _ = f.write(f"The delimited file associating CXX barcodes with samplenames: {samplefile}\n")
            _ = f.write(f"QC checks were performed on demultiplexed FASTQ files using:\n")
            _ = f.write(f"    falco -skip-report -skip-summary input.fq.gz\n")

rule all:
    default_target: True
    input:
        fw_reads = expand(outdir + "{sample}.F.fq.gz", sample = samplenames),
        rv_reads = expand(outdir + "{sample}.R.fq.gz", sample = samplenames),
        runlog   = outdir + "logs/harpy.demultiplex.log",
        qcreport = outdir + "logs/demultiplex.QC.html"
    message:
        "Demultiplexing has finished!"
