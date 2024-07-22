containerized: "docker://pdimens/harpy:latest"

import os
import re
import sys
import multiprocessing
from pathlib import Path
from rich import print as rprint
from rich.panel import Panel

envdir      = os.getcwd() + "/.harpy_envs"
genomefile  = config["inputs"]["genome"]
bamlist     = config["inputs"]["alignments"]
vcffile     = config["inputs"]["vcf"]
samplenames = {Path(i).stem for i in bamlist}
extra       = config.get("extra", None) 
mol_dist    = config["molecule_distance"]
min_sv      = config["min_sv"]
min_barcodes = config["min_barcodes"]
outdir      = config["output_directory"]
skipreports = config["skip_reports"]
bn          = os.path.basename(genomefile)
if bn.lower().endswith(".gz"):
    validgenome = bn[:-3]
else:
    validgenome = bn

if vcffile.lower().endswith("bcf"):
    vcfindex = vcffile + ".csi"
else:
    vcfindex = vcffile + ".tbi"

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

onerror:
    print("")
    rprint(
        Panel(
            f"The workflow has terminated due to an error. See the log file below for more details.",
            title = "[bold]harpy sv naibr",
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
            title = "[bold]harpy sv naibr",
            title_align = "left",
            border_style = "green"
            ),
        file = sys.stderr
    )

def process_args(args):
    argsDict = {
        "min_mapq" : 30,
        "d"        : mol_dist,
        "min_sv"   : min_sv,
        "k"        : min_barcodes
    }
    if args:
        words = [i for i in re.split(r"\s|=", args) if len(i) > 0]
        for i in zip(words[::2], words[1::2]):
            if "blacklist" in i or "candidates" in i:
                argsDict[i[0].lstrip("-")] = i[1]
    return argsDict

def sam_index(infile):
    """Use Samtools to index an input file, adding .bai to the end of the name"""
    if not os.path.exists(f"{infile}.bai"):
        subprocess.run(f"samtools index {infile} {infile}.bai".split())

def get_alignments(wildcards):
    """returns a list with the bam file for the sample based on wildcards.sample"""
    r = re.compile(fr".*/({wildcards.sample})\.(bam|sam)$", flags = re.IGNORECASE)
    aln = list(filter(r.match, bamlist))
    return aln[0]

def get_align_index(wildcards):
    """returns a list with the bai index file for the sample based on wildcards.sample"""
    r = re.compile(fr"(.*/{wildcards.sample})\.(bam|sam)$", flags = re.IGNORECASE)
    aln = list(filter(r.match, bamlist))
    return aln[0] + ".bai"

rule genome_link:
    input:
        genomefile
    output: 
        f"Genome/{validgenome}"
    container:
        None
    message: 
        "Preprocessing {input}"
    shell: 
        """
        if (file {input} | grep -q compressed ) ;then
            # decompress gzipped
            gzip -d -c {input} | seqtk seq > {output}
        elif (file {input} | grep -q BGZF ); then
            # decompress bgzipped
            gzip -d -c {input} | seqtk seq > {output}
        else
            # copy uncompressed
            cp -f {input} {output}
        fi
        """

rule genome_faidx:
    input: 
        f"Genome/{validgenome}"
    output: 
        f"Genome/{validgenome}.fai"
    log:
        f"Genome/{validgenome}.faidx.log"
    container:
        None
    message:
        "Indexing {input}"
    shell:
        "samtools faidx --fai-idx {output} {input} 2> {log}"

rule index_original_alignments:
    input:
        bamlist
    output:
        [f"{i}.bai" for i in bamlist]
    threads:
        workflow.cores
    message:
        "Indexing alignment files"
    run:
        with multiprocessing.Pool(processes=threads) as pool:
            pool.map(sam_index, input)

rule index_bcf:
    input:
        vcffile
    output:
        vcffile + ".csi"
    container:
        None
    message:
        "Indexing {input}"
    shell:
        "bcftools index {input}"

rule index_vcfgz:
    input:
        vcffile
    output:
        vcffile + ".tbi"
    container:
        None
    message:
        "Indexing {input}"
    shell:
        "tabix {input}"

rule phase_alignments:
    input:
        get_align_index,
        vcfindex,
        f"Genome/{validgenome}.fai",
        vcf = vcffile,
        aln = get_alignments,
        ref = f"Genome/{validgenome}"
    output:
        bam = outdir + "/phasedbam/{sample}.bam",
        log = outdir + "/logs/whatshap-haplotag/{sample}.phase.log"
    conda:
        f"{envdir}/phase.yaml"
    threads:
        4
    message:
        "Phasing: {input.aln}"
    shell:
        "whatshap haplotag --sample {wildcards.sample} --ignore-read-groups --tag-supplementary --output-threads={threads} -o {output.bam} --reference {input.ref} {input.vcf} {input.aln} 2> {output.log}"

rule log_phasing:
    input:
        collect(outdir + "/logs/whatshap-haplotag/{sample}.phase.log", sample = samplenames)
    output:
        outdir + "/logs/whatshap-haplotag.log"
    container:
        None
    message:
        "Creating log of alignment phasing"
    shell:
        """
        echo -e "sample\\ttotal_alignments\\tphased_alignments" > {output}
        for i in {input}; do
            SAMP=$(basename $i .phaselog)
            echo -e "${{SAMP}}\\t$(grep "Total alignments" $i)\\t$(grep "could be tagged" $i)" |
                sed 's/ \\+ /\\t/g' | cut -f1,3,5 >> {output}
        done
        """

rule create_config:
    input:
        outdir + "/phasedbam/{sample}.bam"
    output:
        outdir + "/workflow/input/{sample}.naibr"
    params:
        lambda wc: wc.get("sample"),
        min(10, workflow.cores)
    message:
        "Creating naibr config file: {wildcards.sample}"
    run:
        argdict = process_args(extra)
        with open(output[0], "w") as conf:
            _ = conf.write(f"bam_file={input[0]}\n")
            _ = conf.write(f"prefix={params[0]}\n")
            _ = conf.write(f"outdir={outdir}/{params[0]}\n")
            _ = conf.write(f"threads={params[1]}\n")
            for i in argdict:
                _ = conf.write(f"{i}={argdict[i]}\n")

rule index_phased_alignment:
    input:
        outdir + "/phasedbam/{sample}.bam"
    output:
        outdir + "/phasedbam/{sample}.bam.bai"
    container:
        None
    message:
        "Indexing alignment: {wildcards.sample}"
    shell:
        "samtools index {input} {output} 2> /dev/null"

rule call_sv:
    input:
        bam   = outdir + "/phasedbam/{sample}.bam",
        bai   = outdir + "/phasedbam/{sample}.bam.bai",
        conf  = outdir + "/workflow/input/{sample}.naibr"
    output:
        bedpe = outdir + "/{sample}/{sample}.bedpe",
        refmt = outdir + "/{sample}/{sample}.reformat.bedpe",
        vcf   = outdir + "/{sample}/{sample}.vcf"
    log:
        outdir + "/logs/{sample}.naibr.log"
    threads:
        10
    conda:
        f"{envdir}/sv.yaml"
    message:
        "Calling variants: {wildcards.sample}"
    shell:
        "naibr {input.conf} > {log} 2>&1"

rule infer_sv:
    input:
        bedpe = outdir + "/{sample}/{sample}.bedpe",
        refmt = outdir + "/{sample}/{sample}.reformat.bedpe",
        vcf   = outdir + "/{sample}/{sample}.vcf"
    output:
        bedpe = outdir + "/bedpe/{sample}.bedpe",
        refmt = outdir + "/IGV/{sample}.reformat.bedpe",
        fail  = outdir + "/filtered/{sample}.fail.bedpe",
        vcf   = outdir + "/vcf/{sample}.vcf" 
    params:
        outdir = lambda wc: outdir + "/" + wc.get("sample")
    container:
        None
    message:
        "Inferring variants from naibr output: {wildcards.sample}"
    shell:
        """
        infer_sv.py {input.bedpe} -f {output.fail} > {output.bedpe}
        mv {input.refmt} {output.refmt} &&
        mv {input.vcf} {output.vcf} &&
        rm -rf {params.outdir}
        """

rule merge_variants:
    input:
        collect(outdir + "/bedpe/{sample}.bedpe", sample = samplenames)
    output:
        outdir + "/inversions.bedpe",
        outdir + "/deletions.bedpe",
        outdir + "/duplications.bedpe"
    message:
        "Aggregating the detected variants"
    run:
        from pathlib import Path
        with open(output[0], "w") as inversions, open(output[1], "w") as deletions, open(output[2], "w") as duplications:
            header = ["Sample","Chr1","Break1","Chr2","Break2","SplitMolecules","DiscordantReads","Orientation","Haplotype","Score","PassFilter","SV"]
            _ = inversions.write("\t".join(header) + "\n")
            _ = deletions.write("\t".join(header) + "\n")
            _ = duplications.write("\t".join(header) + "\n")
            for varfile in input:
                samplename = Path(varfile).stem
                with open(varfile, "r") as f_in:
                    # read the header to skip it
                    f_in.readline()
                    # read the rest of it
                    while True:
                        line = f_in.readline()
                        if not line:
                            break
                        record = line.rstrip().split("\t")
                        if record[-1] == "inversion":
                            _ = inversions.write(f"{samplename}\t{line}")
                        elif record[-1] == "deletion":
                            _ = deletions.write(f"{samplename}\t{line}")
                        elif record[-1] == "duplication":
                            _ = duplications.write(f"{samplename}\t{line}")

rule create_report:
    input:
        bedpe = outdir + "/{sample}.bedpe",
        fai   = f"Genome/{validgenome}.fai"
    output:
        outdir + "/reports/{sample}.naibr.html"
    conda:
        f"{envdir}/r.yaml"
    message:
        "Creating report: {wildcards.sample}"
    script:
        "report/naibr.Rmd"

rule workflow_summary:
    default_target: True
    input:
        bedpe = collect(outdir + "/bedpe/{sample}.bedpe", sample = samplenames),
        bedpe_agg = collect(outdir + "/{sv}.bedpe", sv = ["inversions", "deletions","duplications"]),
        phaselog = outdir + "/logs/whatshap-haplotag.log",
        reports =  collect(outdir + "/reports/{sample}.naibr.html", sample = samplenames) if not skipreports else []
    message:
        "Summarizing the workflow: {output}"
    run:
        os.system(f"rm -rf {outdir}/naibrlog")
        argdict = process_args(extra)
        with open(outdir + "/workflow/sv.naibr.summary", "w") as f:
            _ = f.write("The harpy sv naibr workflow ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write("The alignment files were phased using:\n")
            _ = f.write(f"    whatshap haplotag --reference genome.fasta --linked-read-distance-cutoff {mol_dist} --ignore-read-groups --tag-supplementary --sample sample_x file.vcf sample_x.bam\n")
            _ = f.write("naibr variant calling ran using these configurations:\n")
            _ = f.write(f"    bam_file=BAMFILE\n")
            _ = f.write(f"    prefix=PREFIX\n")
            _ = f.write(f"    outdir=Variants/naibr/PREFIX\n")
            for i in argdict:
                _ = f.write(f"    {i}={argdict[i]}\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")