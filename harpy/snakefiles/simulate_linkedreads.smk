containerized: "docker://pdimens/harpy:latest"

import os
import logging

output_pref = config["output-prefix"]
in_fasta = config["inputs"]["haplotypes"]
in_bc = config["inputs"]["barcodes"]
haps = [f"{i}".zfill(3) for i in range(1,len(in_fasta) + 1)]

onstart:
    logfile_handler = logger_manager._default_filehandler(config["snakemake_log"])
    logger.addHandler(logfile_handler)

rule simulate_reads:
    input:
        in_fasta,
        in_bc if os.path.exists(in_bc) else [],
    output:
        collect(output_pref + "_S1_L{hap}_R{FR}_001.fq.gz", hap = haps, FR = [1,2])
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
        f'-n {config ["molecule-number"]}',
        f'--mutation {config["mutation"]}',
        f'--indels {config["indels"]}',
        f'--extindels {config["extindels"]}',
        f'-o {config["output-prefix"]}',
        f'-O {config["output-type"]}',
        bc = in_bc if not os.path.exists(in_bc) else ""
    threads:
        workflow.cores
    conda:
        "envs/simulations.yaml"
    shell:
        "mimick {params} {input} 2> {log}"
