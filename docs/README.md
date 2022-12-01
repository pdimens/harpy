![logo](_media/harpy.svg)

Harpy is a haplotagging data processing pipeline for Linux-based systems. It uses all the 
magic of [Snakemake](https://snakemake.readthedocs.io/en/stable/) under the hood to handle 
the worklfow decision-making, but as a user, you just interact with it like a normal command-line 
program. Harpy uses both well known and niche programs to take raw haplotagging sequences and process
them to become called SNP genotypes.

Harpy is modular, meaning you can use different parts of it independent from each other. Need to only align reads?
Great! Only want to call variants? Awesome!

Getting started with Harpy is as simple as installing it from conda/mamba (not implemented yet)
```bash
mamba install -c bioconda -c conda-forge harpy
```

## Usage
You can call `harpy` without any arguments (or with `--help`) to print the docstring to your terminal.
```
harpy --help
                                                           
 Usage: harpy [OPTIONS] COMMAND [ARGS]...                     
                                                              
 HARPY Haplotagging data processing pipeline.                 
 The pipeline trims reads, map sequences, calls variants,     
 imputes genotypes, and phases haplotypes. Get started by     
 running harpy init to generate a configuration file and      
 modify it to your needs. The workflow is:                    
                                                              
 init 🡒 trim 🡒 align 🡒 callvariants 🡒 impute 🡒 phase           
                                                              
 Documentation: https://harpy.github.io                       
                                                              
╭─ Options ──────────────────────────────────────────────────╮
│ --help      Show this message and exit.                    │
╰────────────────────────────────────────────────────────────╯
╭─ Commands ─────────────────────────────────────────────────╮
│ align         Align sample sequences to a reference genome │
│ callvariants  Call variants from sample alignments         │
│ impute        Impute genotypes from genotype likelihoods   │
│ init          Generate template configuration file         │
│ phase         Phase SNPs into haplotypes                   │
│ trim          Remove adapters and quality trim sequences   │
╰────────────────────────────────────────────────────────────╯
```

You can likewise call any of the modules with `--help` (e.g. `harpy align`) to see their usage.
```
harpy align --help

 Usage: harpy align [OPTIONS]                                              
                                                                           
 Align sample sequences to a reference genome                              
 If you don't have a configuration file, use harpy init to generate one    
 and modify it for your project.                                           
                                                                           
╭─ Options ───────────────────────────────────────────────────────────────╮
│ --config   -c  PATH     HARPY configuration yaml file                   │
│                         [default: config.yaml]                          │
│ --dir      -d  PATH     Directory with sample sequences                 │
│                         [default: SeqTrimmed]                           │
│ --threads  -t  INTEGER  Number of threads to use                        │
│                         [default: 4]                                    │
│ --bwa      -b           Use BWA MEM (ignores bardcodes) instead of EMA  │
│ --resume   -r           Resume an incomplete run                        │
│ --help                  Show this message and exit.                     │
╰─────────────────────────────────────────────────────────────────────────╯
```