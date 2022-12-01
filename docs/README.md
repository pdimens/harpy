![logo](_media/harpy.svg)

Harpy is a haplotagging data processing pipeline for Linux-based systems. It uses all the 
magic of [Snakemake](https://snakemake.readthedocs.io/en/stable/) under the hood to handle 
the worklfow decision-making, but as a user, you just interact with it like a normal command-line 
program. Harpy uses both well known and niche programs to take raw haplotagging sequences and process
them to become called SNP genotypes. Most of the settings are pre-configured, and the settings you
can modify can be configured in a pre-generated configuration file. There aren't too many, which should
make things a little simpler. 

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

You can likewise call any of the modules with `--help` (e.g. `harpy align --help`) to see their usage.