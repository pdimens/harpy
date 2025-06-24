containerized: "docker://pdimens/harpy:latest"

import os
import logging

output_pref = config["output-prefix"]
in_fasta = config["inputs"]["haplotypes"]
in_bc = config["inputs"]["barcodes"]
read_params = config["read_params"]
lr_params = config["linked_read_params"]
variant_params = config["variant_params"]
circular = '-C' if lr_params["circular"] else ''
out_type = f'-O {lr_params["output-type"]}' if lr_params.get('output-type', None) else ''
seed = config.get("random_seed", None)
onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake"]["log"])
    logger.addHandler(logfile_handler)

rule simulate_reads:
    input:
        in_fasta
    output:
        collect(output_pref + ".R{FR}.fq.gz", FR = [1,2])
    log:
        "logs/mimick.simulation.log"        
    params:
        f'--coverage {read_params["read_coverage"]}',
        f'--distance {read_params["outer_distance"]}',
        f'--error {read_params["error_rate"]}',
        f'--lengths {read_params["lengths"]}',
        f'--stdev {read_params["stdev"]}',
        f'-x {lr_params["segments"]}',
        f'-a {lr_params["molecule-attempts"]}',
        f'-c {lr_params["molecule-coverage"]}',
        circular,
        f'-m {lr_params["molecule-length"]}',
        f'-n {lr_params ["molecules-per"]}',
        f'-s {lr_params["singletons"]}',
        f'-S {seed}' if seed else "",
        f'--mutation {variant_params["mutation"]}',
        f'--indels {variant_params["indels"]}',
        f'--extindels {variant_params["extindels"]}',
        f'-o {config["output-prefix"]}',
        out_type,
        bc = in_bc
    threads:
        workflow.cores
    conda:
        "envs/simulations.yaml"
    shell:
        "mimick -q 1 {params} {input} 2> {log}"

rule workflow_summary:
    default_target: True
    input:
        collect(output_pref + ".R{FR}.fq.gz", FR = [1,2])
    params:
        f'--coverage {read_params["read_coverage"]}',
        f'--distance {read_params["outer_distance"]}',
        f'--error {read_params["error_rate"]}',
        f'--lengths {read_params["lengths"]}',
        f'--stdev {read_params["stdev"]}',
        f'-x {lr_params["segments"]}',
        f'-a {lr_params["molecule-attempts"]}',
        f'-c {lr_params["molecule-coverage"]}',
        circular,
        f'-m {lr_params["molecule-length"]}',
        f'-n {lr_params ["molecules-per"]}',
        f'-s {lr_params["singletons"]}',
        f'-S {seed}' if seed else "",
        f'--mutation {variant_params["mutation"]}',
        f'--indels {variant_params["indels"]}',
        f'--extindels {variant_params["extindels"]}',
        f'-o {config["output-prefix"]}',
        out_type,
        bc = in_bc
    run:
        summary = ["The harpy simulate linkedreads workflow ran using these parameters:"]
        mimick = "Mimick ran using:\n"
        mimick += f"\tmimick -q 1 {params} inputs.fasta"
        summary.append(mimick)
        sm = "The Snakemake workflow was called via command line:\n"
        sm += f"\t{config['snakemake']['relative']}"
        summary.append(sm)
        with open("workflow/simulate.linkedreads.summary", "w") as f:
            f.write("\n\n".join(summary))
