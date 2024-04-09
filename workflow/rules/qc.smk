import os
import re
import sys
import glob
from rich.panel import Panel
from rich import print as rprint

min_len 	  = config["min_len"]
max_len 	  = config["max_len"]
extra 	  = config.get("extra", "") 
seq_dir   = config["seq_directory"]
skipadapters  = config["adapters"]
skipreports = config["skipreports"]
outdir      = config["output_directory"]
flist = [os.path.basename(i) for i in glob.iglob(f"{seq_dir}/*") if not os.path.isdir(i)]
r = re.compile(r".*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
fqlist = list(filter(r.match, flist))
bn_r = r"[\.\_][RF](?:[12])?(?:\_00[1-9])*\.f(?:ast)?q(?:\.gz)?$"
samplenames = set([re.sub(bn_r, "", i, flags = re.IGNORECASE) for i in fqlist])

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

onsuccess:
    print("")
    rprint(
        Panel(
            f"The workflow has finished successfully! Find the results in [bold]{outdir}/[/bold]",
            title = "[bold]harpy qc",
            title_align = "left",
            border_style = "green"
            ),
        file = sys.stderr
    )

onerror:
    print("")
    rprint(
        Panel(
            f"The workflow has terminated due to an error. See the log file below for more details.",
            title = "[bold]harpy qc",
            title_align = "left",
            border_style = "red"
            ),
        file = sys.stderr
    )

def get_fq1(wildcards):
    # code that returns a list of fastq files for read 1 based on *wildcards.sample* e.g.
    lst = sorted(glob.glob(seq_dir + "/" + wildcards.sample + "*"))
    r = re.compile(r".*[\_\.][FR][1]?(?:\_00[0-9])*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
    fqlist = list(filter(r.match, lst))
    return fqlist

def get_fq2(wildcards):
    # code that returns a list of fastq files for read 2 based on *wildcards.sample*, e.g.
    lst = sorted(glob.glob(seq_dir + "/" + wildcards.sample + "*"))
    r = re.compile(r".*[\_\.][R][2]?(?:\_00[0-9])*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
    fqlist = list(filter(r.match, lst))
    return fqlist

rule qc_fastp:
    input:
        fw   = get_fq1,
        rv   = get_fq2
    output:
        fw   = outdir + "/{sample}.R1.fq.gz",
        rv   = outdir + "/{sample}.R2.fq.gz",
        json = outdir + "/logs/json/{sample}.fastp.json"
    log:
        html = outdir + "/reports/fastp_reports/{sample}.html",
        serr = outdir + "/logs/fastp_logs/{sample}.log"
    params:
        minlen = f"--length_required {min_len}",
        maxlen = f"--max_len1 {max_len}",
        tim_adapters = "--disable_adapter_trimming" if skipadapters else "--detect_adapter_for_pe",
        extra = extra
    threads:
        2
    conda:
        os.getcwd() + "/.harpy_envs/qc.yaml"
    message:
        "Removing adapters + quality trimming: {wildcards.sample}" if not skipadapters else "Quality trimming: {wildcards.sample}" 
    shell: 
        """
        fastp --trim_poly_g --cut_right {params} --thread {threads} -i {input.fw} -I {input.rv} -o {output.fw} -O {output.rv} -h {log.html} -j {output.json} -R "{wildcards.sample} QC Report" 2> {log.serr}
        """

rule count_beadtags:
    input:
        outdir + "/{sample}.R1.fq.gz"
    output: 
        temp(outdir + "/logs/bxcount/{sample}.count.log")
    message:
        "Counting barcode frequency: {wildcards.sample}"
    conda:
        os.getcwd() + "/.harpy_envs/qc.yaml"
    script:
        "scripts/countBX.py"

rule beadtag_counts_summary:
    input: 
        countlog = expand(outdir + "/logs/bxcount/{sample}.count.log", sample = samplenames)
    output:
        outdir + "/reports/barcode.summary.html"
    conda:
        os.getcwd() + "/.harpy_envs/r-env.yaml"
    message:
        "Summarizing sample barcode validation"
    script:
        "report/BxCount.Rmd"
   
rule create_report:
    input: 
        expand(outdir + "/logs/json/{sample}.fastp.json", sample = samplenames)
    output:
        outdir + "/reports/qc.report.html"
    params:
        outdir
    conda:
        os.getcwd() + "/.harpy_envs/qc.yaml"
    message:
        "Aggregating fastp reports"
    shell: 
        """
        multiqc {params}/logs/json -m fastp --force --filename {output} --quiet --title "QC Summary" --comment "This report aggregates trimming and quality control metrics reported by fastp" --no-data-dir 2>/dev/null
        """

rule log_workflow:
    default_target: True
    input:
        fq = expand(outdir + "/{sample}.{FR}.fq.gz", FR = ["R1", "R2"], sample = samplenames),
        bx_report = outdir + "/reports/barcode.summary.html" if not skipreports else [],
        agg_report = outdir + "/reports/qc.report.html" if not skipreports else []
    output:
        output + "/workflow/qc.summary"    
    params:
        minlen = f"--length_required {min_len}",
        maxlen = f"--max_len1 {max_len}",
        tim_adapters = "--disable_adapter_trimming" if skipadapters else "--detect_adapter_for_pe",
        extra = extra
    message:
        "Summarizing the workflow: {output}"
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy qc module ran using these parameters:\n\n")
            _ = f.write(f"The directory with sequences: {seq_dir}\n")
            _ = f.write("fastp trimming ran using:\n")
            _ = f.write("    fastp --trim_poly_g --cut_right " + " ".join(params) + "\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")