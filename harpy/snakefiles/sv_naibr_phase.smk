containerized: "docker://pdimens/harpy:latest"

import os
import re
import logging
from pathlib import Path

onstart:
    logger.logger.addHandler(logging.FileHandler(config["snakemake_log"]))
onsuccess:
    os.remove(logger.logfile)
onerror:
    os.remove(logger.logfile)
wildcard_constraints:
    sample = "[a-zA-Z0-9._-]+"

outdir      = config["output_directory"]
envdir      = os.path.join(os.getcwd(), outdir, "workflow", "envs")
genomefile  = config["inputs"]["genome"]
bamlist     = config["inputs"]["alignments"]
bamdict     = dict(zip(bamlist, bamlist))
vcffile     = config["inputs"]["vcf"]
samplenames = {Path(i).stem for i in bamlist}
extra       = config.get("extra", None) 
mol_dist    = config["molecule_distance"]
min_quality  = config["min_quality"]
min_sv      = config["min_sv"]
min_barcodes = config["min_barcodes"]
plot_contigs = config["reports"]["plot_contigs"]    
skip_reports = config["reports"]["skip"]
bn          = os.path.basename(genomefile)
if bn.lower().endswith(".gz"):
    validgenome = bn[:-3]
else:
    validgenome = bn
if vcffile.lower().endswith("bcf"):
    vcfindex = vcffile + ".csi"
else:
    vcfindex = vcffile + ".tbi"

def process_args(args):
    argsDict = {
        "min_mapq" : min_quality,
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

argdict = process_args(extra)

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

rule process_genome:
    input:
        genomefile
    output: 
        f"Genome/{validgenome}"
    container:
        None
    shell: 
        "seqtk seq {input} > {output}"

rule index_genome:
    input: 
        f"Genome/{validgenome}"
    output: 
        f"Genome/{validgenome}.fai"
    log:
        f"Genome/{validgenome}.faidx.log"
    container:
        None
    shell:
        "samtools faidx --fai-idx {output} {input} 2> {log}"

rule index_alignments:
    input:
        lambda wc: bamdict[wc.bam]
    output:
        "{bam}.bai"
    container:
        None
    shell:
        "samtools index {input}"

rule index_snps:
    input:
        vcffile
    output:
        vcffile + ".csi"
    container:
        None
    shell:
        "bcftools index {input}"

rule index_snps_gz:
    input:
        vcffile
    output:
        vcffile + ".tbi"
    container:
        None
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
    params:
        mol_dist
    conda:
        f"{envdir}/phase.yaml"
    threads:
        4
    shell:
        "whatshap haplotag --sample {wildcards.sample} --linked-read-distance-cutoff {params} --ignore-read-groups --tag-supplementary --output-threads={threads} -o {output.bam} --reference {input.ref} {input.vcf} {input.aln} 2> {output.log}"

rule log_phasing:
    input:
        collect(outdir + "/logs/whatshap-haplotag/{sample}.phase.log", sample = samplenames)
    output:
        outdir + "/logs/whatshap-haplotag.log"
    container:
        None
    shell:
        """
        echo -e "sample\\ttotal_alignments\\tphased_alignments" > {output}
        for i in {input}; do
            SAMP=$(basename $i .phaselog)
            echo -e "${{SAMP}}\\t$(grep "Total alignments" $i)\\t$(grep "could be tagged" $i)" |
                sed 's/ \\+ /\\t/g' | cut -f1,3,5 >> {output}
        done
        """

rule naibr_config:
    input:
        outdir + "/phasedbam/{sample}.bam"
    output:
        outdir + "/workflow/input/{sample}.naibr"
    params:
        lambda wc: wc.get("sample"),
        min(10, workflow.cores - 1)
    run:
        with open(output[0], "w") as conf:
            _ = conf.write(f"bam_file={input[0]}\n")
            _ = conf.write(f"prefix={params[0]}\n")
            _ = conf.write(f"outdir={outdir}/{params[0]}\n")
            _ = conf.write(f"threads={params[1]}\n")
            for i in argdict:
                _ = conf.write(f"{i}={argdict[i]}\n")

rule index_phased:
    input:
        outdir + "/phasedbam/{sample}.bam"
    output:
        outdir + "/phasedbam/{sample}.bam.bai"
    container:
        None
    shell:
        "samtools index {input} {output} 2> /dev/null"

rule call_variants:
    input:
        bam   = outdir + "/phasedbam/{sample}.bam",
        bai   = outdir + "/phasedbam/{sample}.bam.bai",
        conf  = outdir + "/workflow/input/{sample}.naibr"
    output:
        bedpe = temp(outdir + "/{sample}/{sample}.bedpe"),
        refmt = temp(outdir + "/{sample}/{sample}.reformat.bedpe"),
        vcf   = temp(outdir + "/{sample}/{sample}.vcf"),
        log   = temp(outdir + "/{sample}/{sample}.log")
    log:
        outdir + "/logs/naibr/{sample}.naibr.log"
    threads:
        10
    conda:
        f"{envdir}/variants.yaml"
    shell:
        "naibr {input.conf} > {log} 2>&1"

rule infer_variants:
    priority: 100
    input:
        bedpe = outdir + "/{sample}/{sample}.bedpe",
        refmt = outdir + "/{sample}/{sample}.reformat.bedpe",
        vcf   = outdir + "/{sample}/{sample}.vcf"
    output:
        bedpe = outdir + "/bedpe/{sample}.bedpe",
        refmt = outdir + "/IGV/{sample}.reformat.bedpe",
        fail  = outdir + "/bedpe/qc_fail/{sample}.fail.bedpe",
        vcf   = outdir + "/vcf/{sample}.vcf" 
    container:
        None
    shell:
        """
        infer_sv.py {input.bedpe} -f {output.fail} > {output.bedpe}
        cp {input.refmt} {output.refmt}
        cp {input.vcf} {output.vcf}
        """

rule aggregate_variants:
    input:
        collect(outdir + "/bedpe/{sample}.bedpe", sample = samplenames)
    output:
        outdir + "/inversions.bedpe",
        outdir + "/deletions.bedpe",
        outdir + "/duplications.bedpe"
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

rule sample_reports:
    input: 
        faidx = f"Genome/{bn}.fai",
        bedpe = outdir + "/bedpe/{sample}.bedpe",
        qmd   = f"{outdir}/workflow/report/naibr.qmd"
    output:
        report = outdir + "/reports/{sample}.naibr.html",
        qmd = temp(outdir + "/reports/{sample}.naibr.qmd")
    log:
        outdir + "/logs/reports/{sample}.report.log"
    params:
        sample= lambda wc: "-P sample:" + wc.get('sample'),
        contigs= f"-P contigs:{plot_contigs}"
    conda:
        f"{envdir}/r.yaml"
    shell:
        """
        cp {input.qmd} {output.qmd}
        FAIDX=$(realpath {input.faidx})
        BEDPE=$(realpath {input.bedpe})
        quarto render {output.qmd} -P faidx:$FAIDX -P bedpe:$BEDPE {params} 2> {log}
        """

rule workflow_summary:
    default_target: True
    input:
        bedpe = collect(outdir + "/bedpe/{sample}.bedpe", sample = samplenames),
        bedpe_agg = collect(outdir + "/{sv}.bedpe", sv = ["inversions", "deletions","duplications"]),
        phaselog = outdir + "/logs/whatshap-haplotag.log",
        reports =  collect(outdir + "/reports/{sample}.naibr.html", sample = samplenames) if not skip_reports else []
    run:
        os.system(f"rm -rf {outdir}/naibrlog")
        summary = ["The harpy sv naibr workflow ran using these parameters:"]
        summary.append(f"The provided genome: {bn}")
        phase = "The alignment files were phased using:\n"
        phase += f"\twhatshap haplotag --reference genome.fasta --linked-read-distance-cutoff {mol_dist} --ignore-read-groups --tag-supplementary --sample sample_x file.vcf sample_x.bam"
        summary.append(phase)
        naibr = "naibr variant calling ran using these configurations:\n"
        naibr += "\tbam_file=BAMFILE\n"
        naibr += "\tprefix=PREFIX\n"
        naibr += "\toutdir=Variants/naibr/PREFIX\n"
        naibr += "\n\t".join([f"{k}={v}" for k,v in argdict.items()])
        summary.append(naibr)
        sm = "The Snakemake workflow was called via command line:\n"
        sm = f"\t{config['workflow_call']}"
        summary.append(sm)
        with open(outdir + "/workflow/sv.naibr.summary", "w") as f:
            f.write("\n\n".join(summary))