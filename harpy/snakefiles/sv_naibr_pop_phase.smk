containerized: "docker://pdimens/harpy:latest"

import os
import re
import logging
from pathlib import Path

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake"]["log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+",
    population = r"[a-zA-Z0-9._-]+"

genomefile   = config["inputs"]["reference"]
bn           = os.path.basename(genomefile)
bamlist      = config["inputs"]["alignments"]
bamdict      = dict(zip(bamlist, bamlist))
samplenames  = {Path(i).stem for i in bamlist}
groupfile    = config["inputs"]["groupings"]
vcffile      = config["inputs"]["vcf"]
vcfindex     = (vcffile + ".csi") if vcffile.lower().endswith("bcf") else (vcffile + ".tbi")
extra        = config.get("extra", None) 
min_size       = config["min_size"]
min_quality  = config["min_quality"]
min_barcodes = config["min_barcodes"]
mol_dist     = config["molecule_distance"]
skip_reports  = config["reports"]["skip"]
plot_contigs = config["reports"]["plot_contigs"]    
if bn.lower().endswith(".gz"):
    workflow_geno = f"workflow/reference/{bn[:-3]}"
else:
    workflow_geno = f"workflow/reference/{bn}"

def process_args(args):
    argsDict = {
        "min_mapq" : min_quality,
        "d"        : mol_dist,
        "min_sv"   : min_size,
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

rule preprocess_reference:
    input:
        genomefile
    output: 
        geno = workflow_geno,
        fai = f"{workflow_geno}.fai"
    log:
        f"{workflow_geno}.preprocess.log"
    container:
        None
    shell: 
        """
        seqtk seq {input} > {output}
        samtools faidx --fai-idx {output.fai} {output.geno} 2> {log}
        """

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
        f"{workflow_geno}.fai",
        vcf = vcffile,
        aln = get_alignments,
        ref = workflow_geno
    output:
        bam = temp("phasedbam/{sample}.bam"),
        log = "logs/whatshap-haplotag/{sample}.phase.log"
    params:
        mol_dist
    threads:
        4
    conda:
        "envs/phase.yaml"
    shell:
        "whatshap haplotag --sample {wildcards.sample} --linked-read-distance-cutoff {params} --ignore-read-groups --tag-supplementary --output-threads={threads} -o {output.bam} --reference {input.ref} {input.vcf} {input.aln} 2> {output.log}"

rule log_phasing:
    input:
        collect("logs/whatshap-haplotag/{sample}.phase.log", sample = samplenames)
    output:
        "logs/whatshap-haplotag.log"
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

rule preprocess_groups:
    input:
        grp = groupfile
    output:
        grp = "workflow/sample.groups"
    run:
        with open(input.grp, "r") as infile, open(output.grp, "w") as outfile:
            _ = [outfile.write(i) for i in infile.readlines() if not i.lstrip().startswith("#")]

rule concat_list:
    input:
        "workflow/sample.groups"
    output:
        "workflow/pool_samples/{population}.list"
    run:
        with open(output[0], "w") as fout:
            for bamfile in popdict[wildcards.population]:
                _ = fout.write(f"phasedbam/{Path(bamfile).stem}.bam\n")

rule concat_groups:
    input: 
        bamlist  = "workflow/pool_samples/{population}.list",
        bamfiles = lambda wc: collect("phasedbam/{sample}", sample = popdict[wc.population])
    output:
        temp("workflow/input/concat/{population}.unsort.bam")
    log:
        "logs/concat_groups/{population}.concat.log"
    container:
        None
    shell:
        "concatenate_bam.py -b {input.bamlist} > {output} 2> {log}"

rule sort_groups:
    input:
        "workflow/input/concat/{population}.unsort.bam"
    output:
        bam = temp("workflow/input/{population}.bam"),
        bai = temp("workflow/input/{population}.bam.bai")
    log:
        "logs/samtools_sort/{population}.sort.log"
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
        bam = "workflow/input/{population}.bam"
    output:
        cfg = "workflow/config/{population}.naibr"
    params:
        popu = lambda wc: wc.get("population"),
        thd = min(10, workflow.cores - 1)
    run:
        with open(output.cfg, "w") as conf:
            _ = conf.write(f"bam_file={input.bam}\n")
            _ = conf.write(f"outdir={params.popu}\n")
            _ = conf.write(f"prefix={params.popu}\n")
            _ = conf.write(f"threads={params.thd}\n")
            for i in argdict:
                _ = conf.write(f"{i}={argdict[i]}\n")

rule call_variants:
    input:
        bam   = "workflow/input/{population}.bam",
        bai   = "workflow/input/{population}.bam.bai",
        conf  = "workflow/config/{population}.naibr"
    output:
        bedpe = temp("{population}/{population}.bedpe"),
        refmt = temp("{population}/{population}.reformat.bedpe"),
        vcf   = temp("{population}/{population}.vcf"),
        log   = temp("{population}/{population}.log")
    log:
        "logs/naibr/{population}.naibr.log"
    threads:
        min(10, workflow.cores - 1)
    conda:
        "envs/variants.yaml"
    shell:
        "naibr {input.conf} > {log} 2>&1 && rm -rf naibrlog"

rule infer_variants:
    priority: 100
    input:
        bedpe = "{population}/{population}.bedpe",
        refmt = "{population}/{population}.reformat.bedpe",
        vcf   = "{population}/{population}.vcf"
    output:
        bedpe = "bedpe/{population}.bedpe",
        refmt = "IGV/{population}.reformat.bedpe",
        fail  = "bedpe/qc_fail/{population}.fail.bedpe",
        vcf   = "vcf/{population}.vcf" 
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
        collect("bedpe/{population}.bedpe", population = populations)
    output:
        "inversions.bedpe",
        "deletions.bedpe",
        "duplications.bedpe"
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

rule configure_report:
    input:
        yaml = "workflow/report/_quarto.yml",
        scss = "workflow/report/_harpy.scss"
    output:
        yaml = temp("reports/_quarto.yml"),
        scss = temp("reports/_harpy.scss")
    run:
        import shutil
        for i,o in zip(input,output):
            shutil.copy(i,o)

rule sample_reports:
    input: 
        "reports/_quarto.yml",
        "reports/_harpy.scss",
        faidx = f"{workflow_geno}.fai",
        bedpe = "bedpe/{population}.bedpe",
        qmd   = "workflow/report/naibr.qmd"
    output:
        report = "reports/{population}.naibr.html",
        qmd = temp("reports/{population}.naibr.qmd")
    log:
        "logs/reports/{population}.report.log"
    params:
        sample= lambda wc: "-P sample:" + wc.get('population'),
        contigs= f"-P contigs:{plot_contigs}"
    conda:
        "envs/r.yaml"
    retries:
        3
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        FAIDX=$(realpath {input.faidx})
        BEDPE=$(realpath {input.bedpe})
        quarto render {output.qmd} --log {log} --quiet -P faidx:$FAIDX -P bedpe:$BEDPE {params}
        """

rule aggregate_report:
    input: 
        "reports/_quarto.yml",
        "reports/_harpy.scss",
        faidx = f"{workflow_geno}.fai",
        bedpe = collect("bedpe/{pop}.bedpe", pop = populations),
        qmd   = "workflow/report/naibr_pop.qmd"
    output:
        report = "reports/naibr.summary.html",
        qmd = temp("reports/naibr.summary.qmd")
    log:
        "logs/reports/summary.report.log"
    params:
        bedpedir = "bedpe",
        contigs = f"-P contigs:{plot_contigs}"
    conda:
        "envs/r.yaml"
    retries:
        3
    shell:
        """
        cp -f {input.qmd} {output.qmd}
        FAIDX=$(realpath {input.faidx})
        INPATH=$(realpath {params.bedpedir})
        quarto render {output.qmd} --log {log} --quiet -P faidx:$FAIDX -P bedpedir:$INPATH {params.contigs}
        """

rule workflow_summary:
    default_target: True
    input:
        bedpe = collect("bedpe/{pop}.bedpe", pop = populations),
        bedpe_agg = collect("{sv}.bedpe", sv = ["inversions", "deletions","duplications"]),
        phaselog = "logs/whatshap-haplotag.log",
        reports = collect("reports/{pop}.naibr.html", pop = populations) if not skip_reports else [],
        agg_report = "reports/naibr.summary.html" if not skip_reports else []
    run:
        summary = ["The harpy sv naibr workflow ran using these parameters:"]
        summary.append(f"The provided reference genome: {bn}")
        phase = "The alignment files were phased using:\n"
        phase += f"\twhatshap haplotag --reference reference.fasta --linked-read-distance-cutoff {mol_dist} --ignore-read-groups --tag-supplementary --sample sample_x file.vcf sample_x.bam"
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
        sm = f"\t{config['snakemake']['relative']}"
        summary.append(sm)
        with open("workflow/sv.naibr.summary", "w") as f:
            f.write("\n\n".join(summary))
