---
label: Home
description: Using Harpy to process your haplotagged data
icon: home
---

# :icon-home: Home

![](static/logo_trans.png)

Harpy is a haplotagging data processing pipeline for Linux-based systems. It uses all the 
magic of [Snakemake](https://snakemake.readthedocs.io/en/stable/) under the hood to handle 
the worklfow decision-making, but as a user, you just interact with it like a normal command-line 
program. Harpy uses both well known and niche programs to take raw haplotagging sequences and process
them to become called SNP genotypes (or haplotypes) or large structural variants (inversions, deletions, duplications).
Most of the settings are pre-configured and the settings you can modify are done at the command line. 

## What is haplotagging?
Linked-read sequencing exists to combine the throughput and accuracy of short-read
sequencing with the long range haplotype information of long-read sequencing.
Haplotagging is an implementation of linked-read sequencing developed by
[Meier _et al._](https://doi.org/10.1073/pnas.2015005118) to:

1. sequence a large number of samples
2. achieve high molecular resolution
3. do both within a reasonable budget

If you don't have haplotagged data, then Harpy will likely be of little to no use to you. See the [haplotagging site](https://www.fml.tuebingen.mpg.de/9418/haplotagging)
for more information about haplotagging and why you might consider it for your study system.


## Harpy Modules
Harpy is modular, meaning you can use different parts of it independent from each other. Need to only align reads?
Great! Only want to call variants? Awesome! All modules are called by `harpy <module>`. For example, use `harpy align` to align reads.

|      Module                                                        | Description                                   |
|:-------------------------------------------------------------------|:----------------------------------------------|
| [!badge corners="pill" text="preflight"](Modules/preflight.md)     | Run various format checks for FASTQ and BAM files |
| [!badge corners="pill" text="demultiplex"](Modules/demultiplex.md) | Demultiplex haplotagged FASTQ files           |
| [!badge corners="pill" text="qc"](Modules/qc.md)                   | Remove adapters and quality trim sequences    |
| [!badge corners="pill" text="align"](Modules/Align/bwa.md)         | Align sample sequences to a reference genome  |
| [!badge corners="pill" text="snp"](Modules/snp.md)                 | Call SNPs and small indels                    |
| [!badge corners="pill" text="sv"](Modules/SV/naibr.md)             | Call large structural variants                |
| [!badge corners="pill" text="impute"](Modules/impute.md)           | Impute genotypes using variants and sequences |
| [!badge corners="pill" text="phase"](Modules/phase.md)             | Phase SNPs into haplotypes                    |

## Using Harpy
You can call `harpy` without any arguments (or with `--help`) to print the docstring to your terminal. You can likewise call any of the modules without arguments or with `--help` to see their usage  (e.g. `harpy align --help`).
``` harpy --help                                                      
 Usage: harpy COMMAND [ARGS]...                     
                                                              
                  Harpy haplotagging pipeline                  
 An automated workflow to demultiplex sequences, trim and qc   
 reads, map sequences, call variants, impute genotypes, and    
 phase haplotypes of Haplotagging data. Batteries included.    
                                                               
 demultiplex >> qc >> align >> snp >> impute >> phase >> sv          
                                                               
 Documentation: https://pdimens.github.io/harpy/               
                                                               
╭─ Options ──────────────────────────────────────────────────╮
│ --version      Show the version and exit.                  │
│ --help     -h  Show this message and exit.                 │
╰────────────────────────────────────────────────────────────╯
╭─ Modules ──────────────────────────────────────────────────╮
│ demultiplex  Demultiplex haplotagged FASTQ files           │
│ qc           Remove adapters and quality trim sequences    │
│ align        Align sample sequences to a reference genome  │
│ snp          Call SNPs and small indels                    │
│ sv           Call large structural variants                │
│ impute       Impute genotypes using variants and sequences │
│ phase        Phase SNPs into haplotypes                    │
╰────────────────────────────────────────────────────────────╯
╭─ Other Commands ───────────────────────────────────────────╮
│ preflight     Run file format checks on haplotag data      │
│ popgroup      Create a sample grouping file                │
│ stitchparams  Create a template STITCH parameter file      │
│ hpc           Create a config file to run Harpy on an HPC  │
╰────────────────────────────────────────────────────────────╯
```
