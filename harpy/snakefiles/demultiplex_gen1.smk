containerized: "docker://pdimens/harpy:latest"

import os
import logging
import subprocess

outdir = config["output_directory"]
envdir = os.path.join(os.getcwd(), outdir, "workflow", "envs")
samplefile = config["inputs"]["demultiplex_schema"]
skip_reports = config["reports"]["skip"]
keep_unknown = config["keep_unknown"]
n_chunks = min(int(workflow.cores * 2.5), 999)
onstart:
    logger.logger.addHandler(logging.FileHandler(config["snakemake_log"]))
    os.makedirs(f"{outdir}/reports/data", exist_ok = True)
onsuccess:
    os.remove(logger.logfile)
    try:
        os.remove(f"{outdir}/reads.R1.fq.gz")
        os.remove(f"{outdir}/reads.R2.fq.gz")
        os.remove(f"{outdir}/reads.I1.fq.gz")
        os.remove(f"{outdir}/reads.I2.fq.gz")
    except FileNotFoundError:
        pass
onerror:
    os.remove(logger.logfile)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+",
    FR = r"[12]",
    part = r"\d{3}"

def parse_schema(smpl: str, keep_unknown: bool) -> dict:
    d = {}
    with open(smpl, "r") as f:
        for i in f.readlines():
            # a casual way to ignore empty lines or lines with !=2 fields
            try:
                sample, bc = i.split()
                id_segment = bc[0]
                if sample not in d:
                    d[sample] = [bc]
                else:
                    d[sample].append(bc)
            except ValueError:
                continue
    if keep_unknown:
        d["_unknown_sample"] = f"{id_segment}00"
    return d

samples = parse_schema(samplefile, keep_unknown)
samplenames = [i for i in samples]

def setup_chunks(fq1: str, fq2: str, parts: int) -> dict:
    # find the minimum number of reads between R1 and R2 files
    count_r1 = subprocess.check_output(["zgrep", "-c", "-x", '+', fq1])
    count_r2 = subprocess.check_output(["zgrep", "-c", "-x", '+', fq2])
    read_min = min(int(count_r1), int(count_r2))
    chunks_length = read_min // parts
    starts = list(range(1, read_min, max(chunks_length,20)))
    ends = [i-1 for i in starts[1:]]
    # the last end should be -1, which is the "end" in a seqkit range
    ends[-1] = -1
    formatted_parts = [f"{i:03d}" for i in range(1, parts + 1)]
    # format  {part : (start, end)}
    # example {"001": (1, 5000)}
    return dict(zip(formatted_parts, zip(starts, ends)))

chunk_dict = setup_chunks(config["inputs"]["R1"], config["inputs"]["R2"], n_chunks)
fastq_parts = list(chunk_dict.keys())
rule barcode_segments:
    output:
        collect(outdir + "/workflow/segment_{letter}.bc", letter = ["A","C","B","D"])
    params:
        f"{outdir}/workflow"
    container:
        None
    shell:
        "haplotag_acbd.py {params}"

rule link_input:
    input:
        r1 = config["inputs"]["R1"],
        r2 = config["inputs"]["R2"],
        i1 = config["inputs"]["I1"],
        i2 = config["inputs"]["I2"]
    output:
        r1 = temp(f"{outdir}/reads.R1.fq.gz"),
        r2 = temp(f"{outdir}/reads.R2.fq.gz"),
        i1 = temp(f"{outdir}/reads.I1.fq.gz"),
        i2 = temp(f"{outdir}/reads.I2.fq.gz")
    run:
        for i,o in zip(input,output):
            if os.path.exists(o) or os.path.islink(o):
                os.remove(o)
            os.symlink(i, o)

rule partition_reads:
    group: "partition"
    input:
        outdir + "/reads.R{FR}.fq.gz"
    output:
        temp(outdir + "/reads_chunks/reads.R{FR}.part_{part}.fq.gz")
    params:
        lambda wc: f"-r {chunk_dict[wc.part][0]}:{chunk_dict[wc.part][1]}"
    conda:
        f"{envdir}/demultiplex.yaml"
    shell:
        "seqkit range {params} -o {output} {input}"

use rule partition_reads as partition_index with:
    group: "partition"
    input:
        outdir + "/reads.I{FR}.fq.gz"
    output:
        temp(outdir + "/index_chunks/reads.I{FR}.part_{part}.fq.gz")

checkpoint demultiplex:
    group: "partition"
    priority: 100
    input:
        R1 = outdir + "/reads_chunks/reads.R1.part_{part}.fq.gz",
        R2 = outdir + "/reads_chunks/reads.R2.part_{part}.fq.gz",
        I1 = outdir + "/index_chunks/reads.I1.part_{part}.fq.gz",
        I2 = outdir + "/index_chunks/reads.I2.part_{part}.fq.gz",
        segment_a = f"{outdir}/workflow/segment_A.bc",
        segment_b = f"{outdir}/workflow/segment_B.bc",
        segment_c = f"{outdir}/workflow/segment_C.bc",
        segment_d = f"{outdir}/workflow/segment_D.bc",
        schema = samplefile
    output:
        temp(collect(outdir + "/{sample}.{{part}}.R{FR}.fq", sample = samplenames, FR = [1,2])),
        bx_info = temp(f"{outdir}/logs/part.{{part}}.barcodes")
    log:
        f"{outdir}/logs/demultiplex.{{part}}.log"
    params:
        outdir = outdir,
        qxrx = config["include_qx_rx_tags"],
        keep_unknown = keep_unknown,
        part = lambda wc: wc.get("part")
    conda:
        f"{envdir}/demultiplex.yaml"
    script:
        "scripts/demultiplex_gen1.py"

rule merge_partitions:
    input:
        collect(outdir + "/{{sample}}.{part}.R{{FR}}.fq", part = fastq_parts)
    output:
        outdir + "/{sample}.R{FR}.fq.gz"
    log:
        outdir + "/logs/{sample}.{FR}.concat.log"
    container:
        None
    shell:
        "cat {input} | gzip > {output} 2> {log}"

rule merge_barcode_logs:
    input:
        bc = collect(outdir + "/logs/part.{part}.barcodes", part = fastq_parts)
    output:
        log = f"{outdir}/logs/barcodes.log"
    run:
        bc_dict = {}
        for i in input.bc:
            with open(i, "r") as bc_log:
                # skip first row of column names
                _ = bc_log.readline()
                for line in bc_log:
                    barcode,total,correct,corrected = line.split()
                    bc_stats = [int(total), int(correct), int(corrected)]
                    if barcode not in bc_dict:
                        bc_dict[barcode] = bc_stats
                    else:
                        bc_dict[barcode] = list(map(lambda x,y: x+y, bc_stats, bc_dict[barcode]))
        with open(output.log, "w") as f:
            f.write("Barcode\tTotal_Reads\tCorrect_Reads\tCorrected_Reads\n")
            for k,v in bc_dict.items():
                f.write(k + "\t" + "\t".join([str(i) for i in v]) + "\n")

rule assess_quality:
    input:
        outdir + "/{sample}.R{FR}.fq.gz"
    output: 
        outdir + "/reports/data/{sample}.R{FR}.fastqc"
    log:
        outdir + "/logs/{sample}.R{FR}.qc.log"
    threads:
        1
    conda:
        f"{envdir}/qc.yaml"
    shell:
        """
        ( falco --quiet --threads {threads} -skip-report -skip-summary -data-filename {output} {input} ) > {log} 2>&1 ||
cat <<EOF > {output}
##Falco	1.2.4
>>Basic Statistics	fail
#Measure	Value
Filename	{wildcards.sample}.R{wildcards.FR}.fq.gz
File type	Conventional base calls
Encoding	Sanger / Illumina 1.9
Total Sequences	0
Sequences flagged as poor quality	0
Sequence length	0
%GC	0
>>END_MODULE
EOF      
        """

rule report_config:
    output:
        outdir + "/workflow/multiqc.yaml"
    run:
        import yaml
        configs = {
            "sp": {"fastqc/data": {"fn" : "*.fastqc"}},
            "table_sample_merge": {
                "R1": ".R1",
                "R2": ".R2"
            },
            "title": "Quality Assessment of Demultiplexed Samples",
            "subtitle": "This report aggregates the QA results created by falco",
            "report_comment": "Generated as part of the Harpy demultiplex workflow",
            "report_header_info": [
                {"Submit an issue": "https://github.com/pdimens/harpy/issues/new/choose"},
                {"Read the Docs": "https://pdimens.github.io/harpy/"},
                {"Project Homepage": "https://github.com/pdimens/harpy"}
            ]
        }
        with open(output[0], "w", encoding="utf-8") as yml:
            yaml.dump(configs, yml, default_flow_style= False, sort_keys=False, width=float('inf'))

rule quality_report:
    input:
        fqc = collect(outdir + "/reports/data/{sample}.R{FR}.fastqc", sample = samplenames, FR = [1,2]),
        mqc_yaml = outdir + "/workflow/multiqc.yaml"
    output:
        outdir + "/reports/demultiplex.QA.html"
    log:
        f"{outdir}/logs/multiqc.log"
    params:
        options = "--no-version-check --force --quiet --no-data-dir",
        module = " --module fastqc",
        logdir = outdir + "/reports/data/"
    conda:
        f"{envdir}/qc.yaml"
    shell:
        "multiqc --filename {output} --config {input.mqc_yaml} {params} 2> {log}"

rule workflow_summary:
    default_target: True
    input:
        fq = collect(outdir + "/{sample}.R{FR}.fq.gz", sample = samplenames, FR = [1,2]),
        barcode_logs = f"{outdir}/logs/barcodes.log",
        reports = outdir + "/reports/demultiplex.QA.html" if not skip_reports else []
    params:
        R1 = config["inputs"]["R1"],
        R2 = config["inputs"]["R2"],
        I1 = config["inputs"]["I1"],
        I2 = config["inputs"]["I2"]
    run:
        summary = ["The harpy demultiplex workflow ran using these parameters:"]
        summary.append("Linked Read Barcode Design: Generation I")
        inputs = "The multiplexed input files:\n"
        inputs += f"\tread 1: {params.R1}\n"
        inputs += f"\tread 2: {params.R2}\n"
        inputs += f"\tindex 1: {params.I1}\n"
        inputs += f"\tindex 2: {params.I2}"
        inputs += f"Sample demultiplexing schema: {samplefile}"
        summary.append(inputs)
        chunking = "Input data was partitioned into smaller chunks using:\n"
        chunking += "\tseqkit -r start:stop -o output.fq input.fq"
        summary.append(chunking)
        demux = "Samples were demultiplexed using:\n"
        demux += "\tworkflow/scripts/demultiplex_gen1.py"
        summary.append(demux)
        qc = "QC checks were performed on demultiplexed FASTQ files using:\n"
        qc += "\tfalco -skip-report -skip-summary -data-filename output input.fq.gz"
        summary.append(qc)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['workflow_call']}"
        summary.append(sm)
        with open(outdir + "/workflow/demux.gen1.summary", "w") as f:
            f.write("\n\n".join(summary))
