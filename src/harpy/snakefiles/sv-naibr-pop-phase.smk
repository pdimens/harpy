containerized: "docker://pdimens/harpy:latest"

from rich import print as rprint
from rich.panel import Panel
import sys
import os
import re

bam_dir     = config["seq_directory"]
envdir      = os.getcwd() + "/.harpy_envs"
samplenames = config["samplenames"] 
extra       = config.get("extra", None) 
groupfile   = config["groupings"]
genomefile  = config["genomefile"]
molecule_distance = config["molecule_distance"]
min_sv      = config["min_sv"]
min_barcodes = config["min_barcodes"]
outdir      = config["output_directory"]
skipreports = config["skipreports"]

wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+",
    population = "[a-zA-Z0-9._-]+"


bn = os.path.basename(genomefile)
if bn.lower().endswith(".gz"):
    validgenome = bn[:-3]
else:
    validgenome = bn

vcffile = config["vcf"]
if vcffile.lower().endswith("bcf"):
    vcfindex = vcffile + ".csi"
else:
    vcfindex = vcffile + ".tbi"

def process_args(args):
    argsDict = {
        "min_mapq" : 30,
        "d"        : molecule_distance,
        "min_sv"   : min_sv,
        "k"        : min_barcodes
    }
    if args:
        words = [i for i in re.split(r"\s|=", args) if len(i) > 0]
        for i in zip(words[::2], words[1::2]):
            if "blacklist" in i or "candidates" in i:
                argsDict[i[0].lstrip("-")] = i[1]
    return argsDict

# create dictionary of population => filenames
## this makes it easier to set the snakemake rules/wildcards
def pop_manifest(infile, dirn, sampnames):
    d = dict()
    with open(infile) as f:
        for line in f:
            samp, pop = line.rstrip().split()
            if samp.lstrip().startswith("#"):
                continue
            samp = f"{dirn}/phasedbam/{samp}.bam"
            if pop not in d.keys():
                d[pop] = [samp]
            else:
                d[pop].append(samp)
    return d

popdict     = pop_manifest(groupfile, outdir, samplenames)
populations = popdict.keys()

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
            # linked uncompressed
            ln -sr {input} {output}
        fi
        """

rule genome_faidx:
    input: 
        f"Genome/{validgenome}"
    output: 
        f"Genome/{validgenome}.fai"
    container:
        None
    message:
        "Indexing {input}"
    log:
        f"Genome/{validgenome}.faidx.log"
    shell:
        "samtools faidx --fai-idx {output} {input} 2> {log}"

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

rule index_alignments:
    input:
        bam_dir + "/{sample}.bam"
    output:
        bam_dir + "/{sample}.bam.bai"
    container:
        None
    message:
        "Indexing {input}"
    shell:
        "samtools index {input}"

rule phase_alignments:
    input:
        vcfindex,
        bam_dir + "/{sample}.bam.bai",
        f"Genome/{validgenome}.fai",
        vcf = vcffile,
        aln = bam_dir + "/{sample}.bam",
        ref = f"Genome/{validgenome}"
    output:
        bam = outdir + "/phasedbam/{sample}.bam",
        log = outdir + "/logs/whatshap-haplotag/{sample}.phase.log"
    threads:
        4
    conda:
        f"{envdir}/phase.yaml"
    message:
        "Phasing: {input.aln}"
    shell:
        "whatshap haplotag --sample {wildcards.sample} --ignore-read-groups --tag-supplementary --output-threads={threads} -o {output.bam} --reference {input.ref} {input.vcf} {input.aln} 2> {output.log}"

rule log_phasing:
    input:
        collect(outdir + "/logs/whatshap-haplotag/{sample}.phase.log", sample = samplenames)
    output:
        outdir + "/logs/whatshap-haplotag/phasing.log"
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

rule copy_groupings:
    input:
        groupfile
    output:
        outdir + "/workflow/sample.groups"
    message:
        "Logging {input}"
    run:
        with open(input[0], "r") as infile, open(output[0], "w") as outfile:
            _ = [outfile.write(i) for i in infile.readlines() if not i.lstrip().startswith("#")]

rule bam_list:
    input:
        outdir + "/workflow/sample.groups"
    output:
        collect(outdir + "/workflow/{pop}.list", pop = populations)
    message:
        "Creating file lists for each population."
    run:
        for p in populations:
            bamlist = popdict[p]
            with open(f"{outdir}/workflow/{p}.list", "w") as fout:
                for bamfile in bamlist:
                    _ = fout.write(bamfile + "\n")

rule merge_populations:
    input: 
        bamlist  = outdir + "/workflow/{population}.list",
        bamfiles = lambda wc: collect("{sample}", sample = popdict[wc.population]) 
    output:
        bam = temp(outdir + "/workflow/inputpop/{population}.bam"),
        bai = temp(outdir + "/workflow/inputpop/{population}.bam.bai")
    threads:
        2
    container:
        None
    message:
        "Merging alignments: Population {wildcards.population}"
    shell:
        "samtools merge -o {output.bam}##idx##{output.bai} --threads {threads} --write-index -b {input.bamlist}"

rule create_config:
    input:
        outdir + "/workflow/inputpop/{population}.bam"
    output:
        outdir + "/workflow/config/{population}.naibr"
    params:
        lambda wc: wc.get("population"),
        min(10, workflow.cores)
    message:
        "Creating naibr config file: {wildcards.population}"
    run:
        argdict = process_args(extra)
        with open(output[0], "w") as conf:
            _ = conf.write(f"bam_file={input[0]}\n")
            _ = conf.write(f"outdir={outdir}/{params[0]}\n")
            _ = conf.write(f"prefix={params[0]}\n")
            _ = conf.write(f"threads={params[1]}\n")
            for i in argdict:
                _ = conf.write(f"{i}={argdict[i]}\n")

rule call_sv:
    input:
        bam   = outdir + "/workflow/inputpop/{population}.bam",
        bai   = outdir + "/workflow/inputpop/{population}.bam.bai",
        conf  = outdir + "/workflow/config/{population}.naibr"
    output:
        bedpe = outdir + "/{population}/{population}.bedpe",
        refmt = outdir + "/{population}/{population}.reformat.bedpe",
        vcf   = outdir + "/{population}/{population}.vcf"
    log:
        outdir + "/logs/{population}.log"
    threads:
        min(10, workflow.cores)
    conda:
        f"{envdir}/sv.yaml"
    message:
        "Calling variants: {wildcards.population}"
    shell:
        "naibr {input.conf} > {log} 2>&1"

rule infer_sv:
    input:
        bedpe = outdir + "/{population}/{population}.bedpe",
        refmt = outdir + "/{population}/{population}.reformat.bedpe",
        vcf   = outdir + "/{population}/{population}.vcf"
    output:
        bedpe = outdir + "/{population}.bedpe",
        refmt = outdir + "/IGV/{population}.reformat.bedpe",
        fail  = outdir + "/filtered/{population}.fail.bedpe",
        vcf   = outdir + "/vcf/{population}.vcf" 
    params:
        outdir = lambda wc: outdir + "/" + wc.get("population")
    container:
        None
    message:
        "Inferring variants from naibr output: {wildcards.population}"
    shell:
        """
        inferSV.py {input.bedpe} -f {output.fail} > {output.bedpe}
        mv {input.refmt} {output.refmt} &&
        mv {input.vcf} {output.vcf} &&
        rm -rf {params.outdir}
        """

rule create_report:
    input:
        fai   = f"Genome/{validgenome}.fai",
        bedpe = outdir + "/{population}.bedpe"
    output:
        outdir + "/reports/{population}.naibr.html"
    message:
        "Creating report: {wildcards.population}"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/Naibr.Rmd"

rule report_pop:
    input:
        fai   = f"Genome/{validgenome}.fai",
        bedpe = collect(outdir + "/{pop}.bedpe", pop = populations)
    output:
        outdir + "/reports/naibr.pop.summary.html"
    message:
        "Creating summary report"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/NaibrPop.Rmd"

rule log_workflow:
    default_target: True
    input:
        bedpe = collect(outdir + "/{pop}.bedpe", pop = populations),
        reports = collect(outdir + "/reports/{pop}.naibr.html", pop = populations) if not skipreports else [],
        agg_report = outdir + "/reports/naibr.pop.summary.html" if not skipreports else []
    message:
        "Summarizing the workflow: {output}"
    run:
        os.system(f"rm -rf {outdir}/naibrlog")
        argdict = process_args(extra)
        with open(outdir + "/workflow/sv.naibr.summary", "w") as f:
            _ = f.write("The harpy variants sv module ran using these parameters:\n\n")
            _ = f.write(f"The provided genome: {bn}\n")
            _ = f.write(f"The directory with alignments: {bam_dir}\n")
            _ = f.write(f"The sample grouping file: {groupfile}\n\n")
            _ = f.write("The alignment files were phased using:\n")
            _ = f.write(f"    whatshap haplotag --reference genome.fasta --linked-read-distance-cutoff {molecule_distance} --ignore-read-groups --tag-supplementary --sample sample_x file.vcf sample_x.bam\n")
            _ = f.write("naibr variant calling ran using these configurations:\n")
            _ = f.write(f"    bam_file=BAMFILE\n")
            _ = f.write(f"    prefix=PREFIX\n")
            _ = f.write(f"    outdir=Variants/naibr/PREFIX\n")
            for i in argdict:
                _ = f.write(f"    {i}={argdict[i]}\n")
            _ = f.write("\nThe Snakemake workflow was called via command line:\n")
            _ = f.write("    " + str(config["workflow_call"]) + "\n")