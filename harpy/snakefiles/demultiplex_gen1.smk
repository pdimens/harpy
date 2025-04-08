containerized: "docker://pdimens/harpy:latest"

import os
import logging

envdir = os.path.join(os.getcwd(), "workflow", "envs")
samplefile = config["inputs"]["demultiplex_schema"]
skip_reports = config["reports"]["skip"]
keep_unknown = config["keep_unknown"]

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake_log"])
    logger.addHandler(logfile_handler)
wildcard_constraints:
    sample = r"[a-zA-Z0-9._-]+",
    FR = r"[12]",
    part = r"\d{3}"

def parse_schema(smpl, keep_unknown):
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
print(samplenames)
fastq_parts = [f"{i:03d}" for i in range(1, min(workflow.cores, 999) + 1)]

rule barcode_segments:
    output:
        collect("workflow/segment_{letter}.bc", letter = ["A","C","B","D"])
    params:
        "workflow"
    container:
        None
    shell:
        "haplotag_acbd.py {params}"

rule partition_reads:
    input:
        r1 = config["inputs"]["R1"],
        r2 = config["inputs"]["R2"]       
    output:
        r1 = temp("reads.R1.fq.gz"),
        r2 = temp("reads.R2.fq.gz"),
        parts = temp(collect("reads_chunks/reads.R{FR}.part_{part}.fq.gz", part = fastq_parts, FR = [1,2]))
    log:
        "logs/partition.reads.log"
    threads:
        workflow.cores
    params:
        chunks = min(workflow.cores, 999),
        outdir = "reads_chunks"
    conda:
        f"{envdir}/demultiplex.yaml"
    shell:
        """
        ln -sr {input.r1} {output.r1}
        ln -sr {input.r2} {output.r2}
        seqkit split2 -f --quiet -1 {output.r1} -2 {output.r2} -p {params.chunks} -j {threads} -O {params.outdir} -e .gz 2> {log}
        """

use rule partition_reads as partition_index with:
    input:
        r1 = config["inputs"]["I1"],
        r2 = config["inputs"]["I2"]       
    output:
        r1 = temp("reads.I1.fq.gz"),
        r2 = temp("reads.I2.fq.gz"),
        parts = temp(collect("index_chunks/reads.I{FR}.part_{part}.fq.gz", part = fastq_parts, FR = [1,2]))
    log:
        "logs/partition.index.log"
    params:
        chunks = min(workflow.cores, 999),
        outdir = "index_chunks"

rule demultiplex:
    input:
        R1 = "reads_chunks/reads.R1.part_{part}.fq.gz",
        R2 = "reads_chunks/reads.R2.part_{part}.fq.gz",
        I1 = "index_chunks/reads.I1.part_{part}.fq.gz",
        I2 = "index_chunks/reads.I2.part_{part}.fq.gz",
        segment_a = "workflow/segment_A.bc",
        segment_b = "workflow/segment_B.bc",
        segment_c = "workflow/segment_C.bc",
        segment_d = "workflow/segment_D.bc",
        schema = samplefile
    output:
        temp(collect("{sample}.{{part}}.R{FR}.fq", sample = samplenames, FR = [1,2])),
        bx_info = temp("logs/part.{part}.barcodes")
    log:
        "logs/demultiplex.{{part}}.log"
    params:
        outdir = os.getcwd(),
        qxrx = config["include_qx_rx_tags"],
        keep_unknown = keep_unknown,
        part = lambda wc: wc.get("part")
    conda:
        f"{envdir}/demultiplex.yaml"
    script:
        "scripts/demultiplex_gen1.py"

rule merge_partitions:
    input:
        collect("{{sample}}.{part}.R{{FR}}.fq", part = fastq_parts)
    output:
        "{sample}.R{FR}.fq.gz"
    log:
        "logs/{sample}.{FR}.concat.log"
    container:
        None
    shell:
        "cat {input} | gzip > {output} 2> {log}"

rule merge_barcode_logs:
    input:
        bc = collect("logs/part.{part}.barcodes", part = fastq_parts)
    output:
        concat = temp("logs/barcodes.concat"),
        log = "logs/barcodes.log"
    run:
        shell(f"cat {input.bc} | sort -k1,1 > {output.concat}")
        with open(output.concat, "r") as file, open(output.log, "w") as file_out:
            file_out.write("Barcode\tTotal_Reads\tCorrect_Reads\tCorrected_Reads\n")
            prev_bc, prev_1, prev_2, prev_3 = file.readline().split()
            # protect against the headers appearing at the top, just in case
            while prev_bc == "Barcode":
                prev_bc, prev_1, prev_2, prev_3 = file.readline.split()
            for line in file:
                # another redundancy
                if line.startswith("Barcode"):
                    continue
                bc, c1, c2, c3 = line.split()
                if bc != prev_bc:
                    # the barcode is different, write the previous line to file
                    _ = file_out.write(f"{prev_bc}\t{prev_1}\t{prev_2}\t{prev_3}\n")
                    # current becomes the basis for comparison (previous)
                    prev_bc, prev_1, prev_2, prev_3 = bc, c1, c2, c3
                else:
                    # the barcode is the same, sum the count columns
                    # the summed row becomes the basis of comparison
                    prev_bc = bc
                    prev_1 = int(prev_1) + int(c1)
                    prev_2 = int(prev_2) + int(c2)
                    prev_3 = int(prev_3) + int(c3)

rule assess_quality:
    input:
        "{sample}.R{FR}.fq.gz"
    output: 
        "reports/data/{sample}.R{FR}.fastqc"
    log:
        "logs/{sample}.R{FR}.qc.log"
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
        "workflow/multiqc.yaml"
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
        fqc = collect("reports/data/{sample}.R{FR}.fastqc", sample = samplenames, FR = [1,2]),
        mqc_yaml = "workflow/multiqc.yaml"
    output:
        "reports/demultiplex.QA.html"
    log:
        "logs/multiqc.log"
    params:
        options = "--no-version-check --force --quiet --no-data-dir",
        module = " --module fastqc",
        logdir = "reports/data/"
    conda:
        f"{envdir}/qc.yaml"
    shell:
        "multiqc --filename {output} --config {input.mqc_yaml} {params} 2> {log}"

rule workflow_summary:
    default_target: True
    input:
        fq = collect("{sample}.R{FR}.fq.gz", sample = samplenames, FR = [1,2]),
        barcode_logs = "logs/barcodes.log",
        reports = "reports/demultiplex.QA.html" if not skip_reports else []
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
        demux = "Samples were demultiplexed using:\n"
        demux += "\tworkflow/scripts/demultiplex_gen1.py"
        summary.append(demux)
        qc = "QC checks were performed on demultiplexed FASTQ files using:\n"
        qc += "\tfalco -skip-report -skip-summary -data-filename output input.fq.gz"
        summary.append(qc)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['snakemake_command']}"
        summary.append(sm)
        with open("workflow/demux.gen1.summary", "w") as f:
            f.write("\n\n".join(summary))