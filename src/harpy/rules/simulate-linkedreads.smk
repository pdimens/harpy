import os
import sys
import gzip
import shutil
from pathlib import Path

from rich.panel import Panel
from rich import print as rprint

outdir = config["output_directory"]
gen_hap1 = config["genome_hap1"]
gen_hap2 = config["genome_hap2"]

barcodes = config.get("barcodes", None)
if barcodes:
    barcodefile = barcodes
else:
    barcodefile = f"{outdir}/workflow/input/4M-with-alts-february-2016.txt"

onsuccess:
    shutil.rmtree('./_Inline', ignore_errors=True)
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
    shutil.rmtree('./_Inline', ignore_errors=True)
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

rule link_haplotype:
    input:
        gen_hap1
    output: 
        f"{outdir}/workflow/input/hap.0.fasta"
    message: 
        "Decompressing {input}" if gen_hap1.lower().endswith("gz") else "Symlinking {input}"
    run:
        if input[0].lower().endswith("gz"):
            with open(input[0], 'rb') as inf, open(output[0], 'w', encoding='utf8') as outf:
                decom_str = gzip.decompress(inf.read()).decode('utf-8')
                outf.write(decom_str)
        else:
            if not (Path(output[0]).is_symlink() or Path(output[0]).exists()):
                Path(output[0]).symlink_to(Path(input[0]).absolute()) 

use rule link_haplotype as link_second_haplotype with:
    input:
        gen_hap2
    output: 
        f"{outdir}/workflow/input/hap.1.fasta"

rule genome_faidx:
    input:
        outdir + "/workflow/input/hap.{hap}.fasta"
    output: 
        outdir + "/workflow/input/hap.{hap}.fasta.fai"
    message:
        "Indexing haplotype {wildcards.hap}"
    shell:
        "samtools faidx --fai-idx {output} {input}"

if not barcodes:
    rule download_barcodes:
        output:
            barcodefile
        message:
            "Downloading list of standard 10X barcodes"
        run:
            from urllib.request import urlretrieve
            _ = urlretrieve("https://raw.githubusercontent.com/aquaskyline/LRSIM/master/4M-with-alts-february-2016.txt", output[0])

rule create_molecules_hap:
    input:
        outdir + "/workflow/input/hap.{hap}.fasta"
    output:
        temp(multiext(outdir + "/dwgsim.{hap}.12", ".bwa.read1.fastq.gz" ,".bwa.read2.fastq.gz", ".mutations.txt", ".mutations.vcf"))
    log:
        outdir + "/logs/dwgsim.hap.{hap}.log"
    params:
        readpairs = int(config["read_pairs"] * 500000),
        outerdist = config["outer_distance"],
        distsd = config["distance_sd"],
        prefix = lambda wc: outdir + "/dwgsim." + wc.get("hap") + ".12"
    conda:
        os.getcwd() + "/.harpy_envs/simulations.yaml"
    message:
        "Creating reads from {input}"
    shell:
        """
        dwgsim -N {params.readpairs} -e 0.0001,0.0016 -E 0.0001,0.0016 -d {params.outerdist} -s {params.distsd} -1 135 -2 151 -H -y 0 -S 0 -c 0 -R 0 -r 0 -F 0 -o 1 -m /dev/null {input} {params.prefix} 2> {log}
        """

rule interleave_dwgsim_output:
    input:
        collect(outdir + "/dwgsim.{{hap}}.12.bwa.read{rd}.fastq.gz", rd = [1,2]) 
    output:
        outdir + "/dwgsim.{hap}.12.fastq"
    message:
        "Decompressing DWGSIM haplotype {wildcards.hap} output"
    shell:
        "seqtk mergepe {input} > {output}"

rule lrsim:
    default_target: True
    input:
        hap1 = f"{outdir}/dwgsim.0.12.fastq",
        hap2 = f"{outdir}/dwgsim.1.12.fastq",
        fai1 = outdir + "/workflow/input/hap.0.fasta.fai",
        fai2 = outdir + "/workflow/input/hap.1.fasta.fai",
        barcodes = barcodefile
    output:
        collect(outdir + "/sim.{hap}.{ext}", hap = [0,1], ext = ["fp", "manifest"])
    log:
        f"{outdir}/logs/LRSIM.log"
    params:
        lrsim = f"{outdir}/workflow/scripts/LRSIMharpy.pl",
        proj_dir = f"{outdir}",
        outdist  = config["outer_distance"],
        dist_sd  = config["distance_sd"],
        n_pairs  = config["read_pairs"],
        mol_len  = config["molecule_length"],
        parts    = config["partitions"],
        mols_per = config["molecules_per_partition"],
        static = "-o 1 -d 2 -u 4"
    threads:
        workflow.cores
    conda:
        os.getcwd() + "/.harpy_envs/simulations.yaml"
    message:
        "Running LRSIM to generate linked reads from\nhaplotype 1: {input.hap1}\nhaplotype 2: {input.hap2}" 
    shell: 
        """
        perl {params.lrsim} -g {input.hap1},{input.hap2} -p {params.proj_dir}/sim \\
            -b {input.barcodes} -r {params.proj_dir} -i {params.outdist} \\
            -s {params.dist_sd} -x {params.n_pairs} -f {params.mol_len} \\
            -t {params.parts} -m {params.mols_per} -z {threads} {params.static} 2> {log}
        rm -f {params.proj_dir}/.sim.status
        """









#rule all:
#    default_target: True
#    input:
#        collect(outdir + "/dwgsim.{hap}.12.fastq", hap = [0,1])