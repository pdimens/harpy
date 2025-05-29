containerized: "docker://pdimens/harpy:latest"

import os
import logging

output_pref = config["output-prefix"]
in_fasta = config["inputs"]["haplotypes"]
in_bc = config["inputs"]["barcodes"]
haps = [f"{i}".zfill(3) for i in range(1,len(in_fasta) + 1)]

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake"]["log"])
    logger.addHandler(logfile_handler)

rule simulate_reads:
    input:
        in_fasta
    output:
        temp(collect(f"mimick/{output_pref}" + ".hap_{hap}.R{FR}.fq.gz", hap = haps, FR = [1,2]))
    log:
        "logs/mimick.simulation.log"        
    params:
        f'--coverage {config["read_coverage"]}',
        f'--distance {config["outer_distance"]}',
        f'--error {config["error_rate"]}',
        f'--length {config["length"]}',
        f'--stdev {config["stdev"]}',
        f'-l {config["lr-type"]}',
        f'-c {config["molecule-coverage"]}',
        f'-m {config["molecule-length"]}',
        f'-n {config ["molecules-per"]}',
        f'--mutation {config["mutation"]}',
        f'--indels {config["indels"]}',
        f'--extindels {config["extindels"]}',
        f'-o mimick/{config["output-prefix"]}',
        f'-O {config["output-type"]}',
        bc = in_bc
    threads:
        workflow.cores
    conda:
        "envs/simulations.yaml"
    shell:
        "mimick -q 1 {params} {input} 2> {log}"

rule proper_pairing:
    input:
        R1 = f"mimick/{output_pref}" + ".hap_{hap}.R1.fq.gz",
        R2 = f"mimick/{output_pref}" + ".hap_{hap}.R2.fq.gz"
    output:
        R1 = output_pref + ".hap_{hap}.R1.fq.gz",
        R2 = output_pref + ".hap_{hap}.R2.fq.gz"
    log:
        "logs/proper_pair.hap_{hap}.log"
    params:
        "--id-regexp '^(\S+)\/[12]'",
        "--force",
        "-O ."
    conda:
        "envs/simulations.yaml"
    shell:
        "seqkit pair {params} -1 {input.R1} -2 {input.R2} 2> {log}"

rule concatenate_haplotypes:
    input:
        collect(output_pref + ".hap_{hap}.R{{FR}}.fq.gz", hap = haps)
    output:
        output_pref + ".R{FR}.fq.gz"
    container:
        None
    shell:
        "cat {input} > {output}"

rule workflow_summary:
    default_target: True
    input:
        collect(output_pref + ".R{FR}.fq.gz", FR = [1,2]),
    params:
        f'--coverage {config["read_coverage"]}',
        f'--distance {config["outer_distance"]}',
        f'--error {config["error_rate"]}',
        f'--length {config["length"]}',
        f'--stdev {config["stdev"]}',
        f'-l {config["lr-type"]}',
        f'-c {config["molecule-coverage"]}',
        f'-m {config["molecule-length"]}',
        f'-n {config ["molecules-per"]}',
        f'--mutation {config["mutation"]}',
        f'--indels {config["indels"]}',
        f'--extindels {config["extindels"]}',
        f'-o {config["output-prefix"]}',
        f'-O {config["output-type"]}',
        bc = in_bc
    run:
        summary = ["The harpy simulate linkedreads workflow ran using these parameters:"]
        mimick = "Mimick ran using:\n"
        mimick += f"\tmimick -q1 {params} inputs.fasta"
        summary.append(mimick)
        pairing = "Haplotypes were concatenated and proper pairing was enforced with seqkit:"
        pairing += "\tseqkit pair -u --id-regexp '^(\S+)\/[12]' --force -1 reads.R1.fq -2 reads.R2.fq"
        summary.append(pairing)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['snakemake']['relative']}"
        summary.append(sm)
        with open("workflow/simulate.linkedreads.summary", "w") as f:
            f.write("\n\n".join(summary))
