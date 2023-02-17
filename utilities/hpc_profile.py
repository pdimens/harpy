#! /usr/bin/env python3

import argparse
import os

parser = argparse.ArgumentParser(prog = 'HPCProfile',
                    description = 'Write a configuration file for running Snakemake on an HPC cluster')
parser.add_argument('hpc', choices = ['slurm', 'sge'], help = "Type of HPC scheduler")
args = parser.parse_args()

if args.hpc == 'slurm':
    if not os.path.exists('slurm'):
        os.makedirs('slurm')
    with open("slurm/config.yaml", "w") as yml:
        yml.write(
            """cluster:
  mkdir -p logs/{rule} &&
  sbatch
    --partition={resources.partition}
    --qos={resources.qos}
    --cpus-per-task={threads}
    --mem={resources.mem_mb}
    --job-name=smk-{rule}-{wildcards}
    --output=logs/{rule}/{rule}-{wildcards}-%j.out
default-resources:
  - partition=<name-of-default-partition>
  - qos=<name-of-quality-of-service>
  - mem_mb=1000
restart-times: 3
max-jobs-per-second: 10
max-status-checks-per-second: 1
local-cores: 1
latency-wait: 60
jobs: 500
keep-going: True
rerun-incomplete: True
printshellcmds: False
scheduler: greedy
use-conda: False
"""
        )
    print("Created HPC profile \'slurm/config.yaml'. Replace the \'partition\' and \'qos\' placeholders in slurm/config.yaml with options relevant to your system.")

