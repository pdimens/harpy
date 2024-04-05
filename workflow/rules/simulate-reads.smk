import os
import re
import sys
import glob
from rich.panel import Panel
from rich import print as rprint

lrsim_params = "-p " + config["output_dir"] + "/sim"
lrsim_params += " -i " + str(config["outer_distance"])
lrsim_params += " -s " + str(config["distance_sd"])
lrsim_params += " -x " + str(config["read_pairs"])
lrsim_params += " -f " + str(config["molecule_length"])
lrsim_params += " -t " + str(config["partitions"])
lrsim_params += " -m " + str(config["molecules_per_partition"])
outdir = config["output_dir"]
gen_hap1 = config["genome_hap1"]
gen_hap2 = config["genome_hap2"]

#TODO DEAL WITH THIS
barcodes = config.get("barcodes", None)

onsuccess:
    print("")
    rprint(
        Panel(
            f"The workflow has finished successfully! Find the results in [bold]{outdir}/[/bold]. If you want to combine both haplotypes of the forward (or reverse) reads together, you can do so with:\n[blue bold]cat reads_hap{{1..2}}.R1.fq.gz > simulations.R1.fq.gz[/blue bold]",
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

rule prepare_directory:
    input:
        hap1 = gen_hap1,
        hap2 = gen_hap2,
    output:
        hap1 = f"{outdir}/sim.hap.0.clean.fasta",
        hap2 = f"{outdir}/sim.hap.1.clean.fasta" 
    message:
        "Fooling LRSIM into a false sense of security"
    shell:
        """
        mkdir -p dwgsim
        ln -s {input.hap1} {output.hap1}
        ln -s {input.hap2} {output.hap2}
        """

rule genome_faidx:
    input:
        outdir + "/sim.hap.{hap}.clean.fasta"
    output: 
        outdir + "/sim.hap.{hap}.clean.fasta.fai"
    log:
        outdir + "/logs/.{hap}.clean.fasta"
    message:
        "Indexing {input}"
    shell:
        "samtools faidx --fai-idx {output} {input} 2> {log}"


rule lrsim:
    input:
        hap1 = f"{outdir}/sim.hap.0.clean.fasta",
        hap2 = f"{outdir}/sim.hap.1.clean.fasta",
        expand(outdir + "/sim.hap.{hap}.clean.fasta.fai", hap = [0,1]),
        barcodes = barcodefile
    output:
        expand(outdir + "/sim_S1_L00{hap}_R{fr}_001.fastq.gz", hap = [0,1], fr = [0,1]),
        expand(outdir + "/sim.{hap}.{ext}", hap = [0,1], ext = ["fp", "manifest", "sort.manifest"]),
        temp(expand(outdir + "/sim.dwgsim.{hap}.12.fastq", hap = [0,1]))
    log:
        f"{outdir}/logs/LRSIM.log"
    params:
        lrsim = f"{outdir}/workflow/scripts/LRSIM.pl",
        runoptions = lrsim_params
    threads:
        workflow.cores
    conda:
        os.getcwd() + "/.harpy_envs/simulations.yaml"
    message:
        f"Running LRSIM to generate linked reads from\nhaplotype 1: {input.hap1}\nhaplotype 2: {input.hap2}" 
    shell: 
        "{params.lrsim} -g {input.hap1},{input.hap2} {params.runoptions} -z {threads} -o -u 3 > {log}"

rule convert_haplotag:
    input:
        fw = outdir + "/sim_S1_L00{hap}_R1_001.fastq.gz",
        fw = outdir + "/sim_S1_L00{hap}_R2_001.fastq.gz",
        barcodes = "BARCODE FILE"
    output:
        fw = "hap{hap}_haplotag.R1.fq.gz",
        rv = "hap{hap}_haplotag.R2.fq.gz"
    log:
        conversions = outdir + "/workflow/10XtoHaplotag_{num}.txt" 
    message:
        "Converting 10X barcodes to haplotag format"
    script:
        "10xtoBX.py"

rule log_runtime:
    output:
        outdir + "/workflow/simulate.reads.workflow.summary"
    params:
        minlen = f"--length_required {min_len}",
        maxlen = f"--max_len1 {max_len}",
        tim_adapters = "--disable_adapter_trimming" if skipadapters else "--detect_adapter_for_pe",
        extra = extra
    message:
        "Creating record of relevant runtime parameters: {output}"
    run:
        with open(output[0], "w") as f:
            _ = f.write("The harpy simulate reads module ran using these parameters:\n\n")
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