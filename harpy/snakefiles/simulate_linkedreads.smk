containerized: "docker://pdimens/harpy:latest"

import os
import logging

output_pref = config["output-prefix"]
in_fasta = config["inputs"]["haplotypes"]
in_bc = config["inputs"]["barcodes"]
read_params = config["read_params"]
lr_params = config["linked_read_params"]
variant_params = config["variant_params"]
seed = config.get("random_seed", None)
onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake"]["log"])
    logger.addHandler(logfile_handler)

rule simulate_reads:
    input:
        in_fasta
    output:
        temp(collect(f"mimick/{output_pref}" + ".R{FR}.fq.gz", FR = [1,2]))
    log:
        "logs/mimick.simulation.log"        
    params:
        f'--coverage {read_params["read_coverage"]}',
        f'--distance {read_params["outer_distance"]}',
        f'--error {read_params["error_rate"]}',
        f'--length {read_params["length"]}',
        f'--stdev {read_params["stdev"]}',
        f'-l {lr_params["lr-type"]}',
        f'-c {lr_params["molecule-coverage"]}',
        f'-m {lr_params["molecule-length"]}',
        f'-n {lr_params ["molecules-per"]}',
        f'-s {lr_params["singletons"]}',
        f'-S {seed}' if seed else "",
        f'--mutation {variant_params["mutation"]}',
        f'--indels {variant_params["indels"]}',
        f'--extindels {variant_params["extindels"]}',
        f'-o mimick/{config["output-prefix"]}',
        f'-O {lr_params["output-type"]}',
        bc = in_bc
    threads:
        workflow.cores
    conda:
        "envs/simulations.yaml"
    shell:
        "mimick -q 1 {params} {input} 2> {log}"

rule proper_pairing:
    input:
        R1 = f"mimick/{output_pref}" + ".R1.fq.gz",
        R2 = f"mimick/{output_pref}" + ".R2.fq.gz"
    output:
        R1 = output_pref + ".R1.fq.gz",
        R2 = output_pref + ".R2.fq.gz"
    log:
        "logs/proper_pair.log"
    params:
        "--id-regexp '^(\S+)\/[12]'",
        "--force",
        "-O ."
    conda:
        "envs/simulations.yaml"
    shell:
        "seqkit pair {params} -1 {input.R1} -2 {input.R2} 2> {log}"

rule workflow_summary:
    default_target: True
    input:
        collect(output_pref + ".R{FR}.fq.gz", FR = [1,2]),
    params:
        f'--coverage {read_params["read_coverage"]}',
        f'--distance {read_params["outer_distance"]}',
        f'--error {read_params["error_rate"]}',
        f'--length {read_params["length"]}',
        f'--stdev {read_params["stdev"]}',
        f'-l {lr_params["lr-type"]}',
        f'-c {lr_params["molecule-coverage"]}',
        f'-m {lr_params["molecule-length"]}',
        f'-n {lr_params ["molecules-per"]}',
        f'-s {lr_params["singletons"]}',
        f'-S {seed}' if seed else "",
        f'--mutation {variant_params["mutation"]}',
        f'--indels {variant_params["indels"]}',
        f'--extindels {variant_params["extindels"]}',
        f'-o mimick/{config["output-prefix"]}',
        f'-O {lr_params["output-type"]}',
        bc = in_bc
    run:
        summary = ["The harpy simulate linkedreads workflow ran using these parameters:"]
        mimick = "Mimick ran using:\n"
        mimick += f"\tmimick -q 1 {params} inputs.fasta"
        summary.append(mimick)
        pairing = "Haplotypes were concatenated and proper pairing was enforced with seqkit:"
        pairing += "\tseqkit pair -u --id-regexp '^(\S+)\/[12]' --force -1 reads.R1.fq -2 reads.R2.fq"
        summary.append(pairing)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['snakemake']['relative']}"
        summary.append(sm)
        with open("workflow/simulate.linkedreads.summary", "w") as f:
            f.write("\n\n".join(summary))
