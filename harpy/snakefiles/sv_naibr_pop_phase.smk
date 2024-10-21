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
    sample = "[a-zA-Z0-9._-]+",
    population = "[a-zA-Z0-9._-]+"

envdir       = os.path.join(os.getcwd(), ".harpy_envs")
genomefile   = config["inputs"]["genome"]
bn           = os.path.basename(genomefile)
bamlist      = config["inputs"]["alignments"]
bamdict      = dict(zip(bamlist, bamlist))
samplenames  = {Path(i).stem for i in bamlist}
groupfile    = config["inputs"]["groupings"]
vcffile      = config["inputs"]["vcf"]
vcfindex     = (vcffile + ".csi") if vcffile.lower().endswith("bcf") else (vcffile + ".tbi")
extra        = config.get("extra", None) 
min_sv       = config["min_sv"]
min_quality  = config["min_quality"]
min_barcodes = config["min_barcodes"]
mol_dist     = config["molecule_distance"]
outdir       = config["output_directory"]
skip_reports  = config["skip_reports"]
if bn.lower().endswith(".gz"):
    bn = bn[:-3]

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

# create dictionary of population => filenames
## this makes it easier to set the snakemake rules/wildcards
def pop_manifest(groupingfile, filelist):
    d = {}
    with open(groupingfile) as f:
        for line in f:
            samp, pop = line.rstrip().split()
            if samp.lstrip().startswith("#"):
                continue
            r = re.compile(fr".*/({samp.lstrip()})\.(bam|sam)$", flags = re.IGNORECASE)
            sampl = list(filter(r.match, filelist))[0]
            sampl = os.path.basename(sampl)
            #f"{out_dir}/phasedbam/{Path(sampl).stem}.bam"
            if pop not in d.keys():
                d[pop] = [sampl]
            else:
                d[pop].append(sampl)
    return d

popdict     = pop_manifest(groupfile, bamlist)
populations = popdict.keys()

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
        f"Genome/{bn}"
    container:
        None
    shell: 
        "seqtk seq {input} > {output}"

rule index_genome:
    input: 
        f"Genome/{bn}"
    output: 
        f"Genome/{bn}.fai"
    container:
        None
    log:
        f"Genome/{bn}.faidx.log"
    shell:
        "samtools faidx --fai-idx {output} {input} 2> {log}"

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

rule index_alignments:
    input:
        lambda wc: bamdict[wc.bam]
    output:
        "{bam}.bai"
    container:
        None
    shell:
        "samtools index {input}"

rule phase_alignments:
    input:
        vcfindex,
        get_align_index,
        f"Genome/{bn}.fai",
        vcf = vcffile,
        aln = get_alignments,
        ref = f"Genome/{bn}"
    output:
        bam = temp(outdir + "/phasedbam/{sample}.bam"),
        log = outdir + "/logs/whatshap-haplotag/{sample}.phase.log"
    params:
        mol_dist
    threads:
        4
    conda:
        f"{envdir}/phase.yaml"
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

rule preproc_groups:
    input:
        groupfile
    output:
        outdir + "/workflow/sample.groups"
    run:
        with open(input[0], "r") as infile, open(output[0], "w") as outfile:
            _ = [outfile.write(i) for i in infile.readlines() if not i.lstrip().startswith("#")]

rule concat_list:
    input:
        outdir + "/workflow/sample.groups"
    output:
        outdir + "/workflow/pool_samples/{population}.list"
    run:
        with open(output[0], "w") as fout:
            for bamfile in popdict[wildcards.population]:
                _ = fout.write(f"{outdir}/phasedbam/{Path(bamfile).stem}.bam\n")

rule concat_groups:
    input: 
        bamlist  = outdir + "/workflow/pool_samples/{population}.list",
        bamfiles = lambda wc: collect(outdir + "/phasedbam/{sample}", sample = popdict[wc.population])
    output:
        temp(outdir + "/workflow/input/concat/{population}.unsort.bam")
    log:
        outdir + "/logs/concat_groups/{population}.concat.log"
    container:
        None
    shell:
        "concatenate_bam.py -o {output} -b {input.bamlist} 2> {log}"

rule sort_groups:
    input:
        outdir + "/workflow/input/concat/{population}.unsort.bam"
    output:
        bam = temp(outdir + "/workflow/input/{population}.bam"),
        bai = temp(outdir + "/workflow/input/{population}.bam.bai")
    log:
        outdir + "/logs/samtools_sort/{population}.sort.log"
    resources:
        mem_mb = 2000
    threads:
        10
    container:
        None
    shell:
        "samtools sort -@ {threads} -O bam -l 0 -m {resources.mem_mb}M --write-index -o {output.bam}##idx##{output.bai} {input} 2> {log}"

rule naibr_config:
    input:
        outdir + "/workflow/input/{population}.bam"
    output:
        outdir + "/workflow/config/{population}.naibr"
    params:
        lambda wc: wc.get("population"),
        min(10, workflow.cores - 1)
    run:
        with open(output[0], "w") as conf:
            _ = conf.write(f"bam_file={input[0]}\n")
            _ = conf.write(f"outdir={outdir}/{params[0]}\n")
            _ = conf.write(f"prefix={params[0]}\n")
            _ = conf.write(f"threads={params[1]}\n")
            for i in argdict:
                _ = conf.write(f"{i}={argdict[i]}\n")

rule call_variants:
    input:
        bam   = outdir + "/workflow/input/{population}.bam",
        bai   = outdir + "/workflow/input/{population}.bam.bai",
        conf  = outdir + "/workflow/config/{population}.naibr"
    output:
        bedpe = outdir + "/{population}/{population}.bedpe",
        refmt = outdir + "/{population}/{population}.reformat.bedpe",
        vcf   = outdir + "/{population}/{population}.vcf"
    log:
        outdir + "/logs/naibr/{population}.naibr.log"
    threads:
        min(10, workflow.cores - 1)
    conda:
        f"{envdir}/variants.yaml"
    shell:
        "naibr {input.conf} > {log} 2>&1"

rule infer_variants:
    input:
        bedpe = outdir + "/{population}/{population}.bedpe",
        refmt = outdir + "/{population}/{population}.reformat.bedpe",
        vcf   = outdir + "/{population}/{population}.vcf"
    output:
        bedpe = outdir + "/bedpe/{population}.bedpe",
        refmt = outdir + "/IGV/{population}.reformat.bedpe",
        fail  = outdir + "/bedpe/qc_fail/{population}.fail.bedpe",
        vcf   = outdir + "/vcf/{population}.vcf" 
    params:
        outdir = lambda wc: outdir + "/" + wc.get("population")
    container:
        None
    shell:
        """
        infer_sv.py {input.bedpe} -f {output.fail} > {output.bedpe}
        mv {input.refmt} {output.refmt} &&
        mv {input.vcf} {output.vcf} &&
        rm -rf {params.outdir}
        """

rule aggregate_variants:
    input:
        collect(outdir + "/bedpe/{population}.bedpe", population = populations)
    output:
        outdir + "/inversions.bedpe",
        outdir + "/deletions.bedpe",
        outdir + "/duplications.bedpe"
    run:
        from pathlib import Path
        with open(output[0], "w") as inversions, open(output[1], "w") as deletions, open(output[2], "w") as duplications:
            header = ["Population","Chr1","Break1","Chr2","Break2","SplitMolecules","DiscordantReads","Orientation","Haplotype","Score","PassFilter","SV"]
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

rule group_reports:
    input:
        fai   = f"Genome/{bn}.fai",
        bedpe = outdir + "/bedpe/{population}.bedpe"
    output:
        outdir + "/reports/{population}.naibr.html"
    log:
        logfile = outdir + "/logs/reports/{population}.report.log"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/naibr.Rmd"

rule aggregate_report:
    input:
        fai   = f"Genome/{bn}.fai",
        bedpe = collect(outdir + "/bedpe/{pop}.bedpe", pop = populations)
    output:
        outdir + "/reports/naibr.pop.summary.html"
    log:
        logfile = outdir + "/logs/reports/summary.report.log"
    conda:
        f"{envdir}/r.yaml"
    script:
        "report/naibr_pop.Rmd"

rule workflow_summary:
    default_target: True
    input:
        bedpe = collect(outdir + "/bedpe/{pop}.bedpe", pop = populations),
        bedpe_agg = collect(outdir + "/{sv}.bedpe", sv = ["inversions", "deletions","duplications"]),
        phaselog = outdir + "/logs/whatshap-haplotag.log",
        reports = collect(outdir + "/reports/{pop}.naibr.html", pop = populations) if not skip_reports else [],
        agg_report = outdir + "/reports/naibr.pop.summary.html" if not skip_reports else []
    run:
        os.system(f"rm -rf {outdir}/naibrlog")
        summary = ["The harpy sv naibr workflow ran using these parameters:"]
        summary.append(f"The provided genome: {bn}")
        phase = "The alignment files were phased using:\n"
        phase += f"\twhatshap haplotag --reference genome.fasta --linked-read-distance-cutoff {mol_dist} --ignore-read-groups --tag-supplementary --sample sample_x file.vcf sample_x.bam"
        summary.append(phase)
        concat = "The alignments were concatenated using:\n"
        concat += "\tconcatenate_bam.py -o groupname.bam -b samples.list"
        summary.append(concat)
        naibr = "naibr variant calling ran using these configurations:\n"
        naibr += "\tbam_file=BAMFILE\n"
        naibr += "\tprefix=PREFIX\n"
        naibr += "\toutdir=Variants/naibr/PREFIX\n"
        naibr += "\n\t".join([f"{k}={v}" for k,v in argdict.items()])
        summary.append(naibr)
        sm = "The Snakemake workflow was called via command line:\n"
        sm = f"\t{config["workflow_call"]}"
        summary.append(sm)
        with open(outdir + "/workflow/sv.naibr.summary", "w") as f:
            f.write("\n\n".join(summary))
