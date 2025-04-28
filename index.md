---
label: Home
description: Using Harpy to process your linked-read data
icon: home
---

# :icon-home: Home

![](static/logo_trans.png)

Harpy is a [haplotagging data](haplotagdata.md) processing pipeline for Linux-based systems-- at least it
was prior to the release of version 2. Now, it can process linked-read data from haplotagging, TELLseq, stLFR, and
even regular non-linked WGS data. It uses all the magic of [Snakemake](https://snakemake.readthedocs.io/en/stable/)
under the hood to handle  the worklfow decision-making, but as a user, you just interact with it like a normal command-line 
program. Harpy employs both well known and niche programs to take raw linked-read sequences and process
them to become called SNP genotypes (or haplotypes) or large structural variants (inversions, deletions, duplications).
Most of the settings are pre-configured and the settings you can modify are done at the command line. Some parts of this documentation
will refer to haplotagging specifically as we either forgot to update parts of the documentation or require you (the user)
to do a data conversion for some parts of Harpy to work with non-haplotagging linked-read data. As always, feel free to drop
an [Issue](https://github.com/pdimens/harpy/issues/new/choose) or open a [Discussion](https://github.com/pdimens/harpy/discussions) on GitHub.

## Harpy Modules
Harpy is modular, meaning you can use different parts of it independent from each other. Need to only align reads?
Great! Only want to call variants? Awesome! All modules are called by `harpy <workflow>`. For example, use `harpy align` to align reads.

{.compact}
| Workflow                                                                | Description                                                          |
| :---------------------------------------------------------------------- | :------------------------------------------------------------------- |
| [!badge corners="pill" text="align"](Workflows/Align/Align.md)          | Align sample sequences to a reference genome                         |
| [!badge corners="pill" text="assembly"](Workflows/assembly.md)          | Create a genome assembly from linked-reads                           |
| [!badge corners="pill" text="demultiplex"](Workflows/demultiplex.md)    | Demultiplex haplotagged FASTQ files                                  |
| [!badge corners="pill" text="impute"](Workflows/impute.md)              | Impute genotypes using variants and sequences                        |
| [!badge corners="pill" text="metassembly"](Workflows/metassembly.md)    | Create a metagenome assembly from linked-reads                       |
| [!badge corners="pill" text="phase"](Workflows/phase.md)                | Phase SNPs into haplotypes                                           |
| [!badge corners="pill" text="preflight"](Workflows/validate.md)         | Run various format checks for FASTQ and BAM files                    |
| [!badge corners="pill" text="qc"](Workflows/qc.md)                      | Remove adapters, deduplicate, and quality trim sequences             |
| [!badge corners="pill" text="simulate"](Workflows/Simulate/Simulate.md) | Simulate haplotag linked reads or genomic variants                   |
| [!badge corners="pill" text="snp"](Workflows/snp.md)                    | Call SNPs and small indels                                           |
| [!badge corners="pill" text="sv"](Workflows/SV/SV.md)                   | Call large structural variants (inversions, deletions, duplications) |

## Using Harpy
You can call `harpy` without any arguments (or with `--help`) to print the docstring to your terminal. You can likewise call any of the modules without arguments or with `--help` to see their usage  (e.g. `harpy align --help`).
``` harpy --help                                                      
 Usage: harpy COMMAND [ARGS]...                                            
                                                                
 An automated workflow for linked-read data to go  
 from raw data to genotypes (or phased haplotypes). Batteries   
 included.                                                      
 demultiplex >> qc >> align >> snp >> impute >> phase >> sv     
                                                                
 Documentation: https://pdimens.github.io/harpy/                
                                                                
╭─ Data Processing ────────────────────────────────────────────╮
│ align        Align sample sequences to a reference genome    │
│ assembly     Create an assembly from linked-reads            │
│ demultiplex  Demultiplex haplotagged FASTQ files             │
│ impute       Impute genotypes using variants and alignments  │
│ metassembly  Create a metassembly from linked-reads          │
│ phase        Phase SNPs into haplotypes                      │
│ qc           Remove adapters and quality-control sequences   │
│ simulate     Simulate variants or linked-reads from a genome │
│ snp          Call SNPs and small indels on alignments        │
│ sv           Call large structural variants on alignments    │
╰──────────────────────────────────────────────────────────────╯
╭─ Other Commands ─────────────────────────────────────────────╮
│ convert     Convert data between linked-read types           │
│ deconvolve  Resolve barcode sharing in unrelated molecules   │
│ downsample  Downsample data by barcode                       │
│ template    Create template files and HPC configs for        │
│             workflows                                        │
╰──────────────────────────────────────────────────────────────╯
╭─ Troubleshoot ───────────────────────────────────────────────╮
│ diagnose  Run the Snakemake debugger to identify hang-ups    │
│ resume    Resume a Harpy workflow from an existing directory │
│ validate  File format checks for linked-read data            │
│ view      View a workflow's components                       │
╰──────────────────────────────────────────────────────────────╯
```

## Typical Linked-Read Workflows
Depending on your project goals, you may want any combination of SNPs, structural
variants (inversions, deletions, duplications), or phased haplotypes. Below is a flow chart
outlining a general workflow of linked-read data.

```mermaid
graph LR
    Demux([demultiplex]):::clean--->QC([QC, trim adapters, etc.]):::clean
    QC--->Align([align sequences]):::clean
    Align--->SNP([call SNPs]):::clean
    SNP--->Impute([impute genotypes]):::clean
    SNP--->Phase([phase haplotypes]):::clean
    Align--->SV([call structural variants]):::clean

    classDef clean fill:#f5f6f9,stroke:#b7c9ef,stroke-width:2px
```

Alternatively, if you are interested in assembling a genome or metagenome, your workflow might look like:

```mermaid
graph LR
    QC([QC, trim adapters, etc.]):::clean--->DC([barcode deconvolution]):::clean
    DC--->Assembly([assembly/metassembly]):::clean

    classDef clean fill:#f5f6f9,stroke:#b7c9ef,stroke-width:2px
```