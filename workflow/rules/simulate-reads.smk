import os
import re
import sys
import glob
from rich.panel import Panel
from rich import print as rprint

onsuccess:
    print("")
    rprint(
        Panel(
            f"The workflow has finished successfully! Find the results in [bold]{outdir}/[/bold]",
            title = "[bold]harpy simulate reads",
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
            title = "[bold]harpy simulate reads",
            title_align = "left",
            border_style = "red"
            ),
        file = sys.stderr
    )

rule lrsim:
    input:
        hap1 = ,
        hap2 = ,
        barcodes = 
    output:
        fw   = outdir + "/{sample}.R1.fq.gz",
        rv   = outdir + "/{sample}.R2.fq.gz"
    log:
        f"{outdir}/logs/LRSIM.log"
    params:
        prefix = outprefix
    threads:
        workflow.cores
    conda:
        os.getcwd() + "/.harpy_envs/simulations.yaml"
    message:
        f"Running LRSIM to generate linked reads from:\nhaplotype 1: {hap1}\nhaplotype 2: {hap2}" 
    shell: 
        """
        
        """

rule log_runtime:
    output:
        outdir + "/workflow/qc.workflow.summary"
    params:
        minlen = f"--length_required {min_len}",
        maxlen = f"--max_len1 {max_len}",
        tim_adapters = "--disable_adapter_trimming" if skipadapters else "--detect_adapter_for_pe",
        extra = extra
    message:
        "Creating record of relevant runtime parameters: {output}"
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy qc module ran using these parameters:\n\n")
            _ = f.write(f"The directory with sequences: {seq_dir}\n")
            _ = f.write("fastp trimming ran using:\n")
            _ = f.write("    fastp --trim_poly_g --cut_right " + " ".join(params) + "\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")

results = list()
results.append(expand(outdir + "/logs/json/{sample}.fastp.json", sample = samplenames))
results.append(expand(outdir + "/{sample}.{FR}.fq.gz", FR = ["R1", "R2"], sample = samplenames))
results.append(outdir + "/workflow/qc.workflow.summary")
if not skipreports:
    results.append(outdir + "/reports/barcode.summary.html")
    
rule create_report:
    default_target: True
    input: 
        results
    output:
        outdir + "/reports/qc.report.html"
    params:
        outdir
    conda:
        os.getcwd() + "/.harpy_envs/qc.yaml"
    message:
        "Sequencing quality filtering and trimming is complete!"
    shell: 
        """
        multiqc {params}/logs/json -m fastp --force --filename {output} --quiet --title "QC Summary" --comment "This report aggregates trimming and quality control metrics reported by fastp" --no-data-dir 2>/dev/null
        """