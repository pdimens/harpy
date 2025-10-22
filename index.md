---
label: Home
description: Using Harpy to process your linked-read data
icon: home
hidden: true
---

# :icon-home: Home

![](static/logo_trans.png)

Harpy is a [haplotagging data](Getting_Started/inputformat.md) processing pipeline for Linux-based systems-- at least it
was prior to the release of version 2. Now, it can process linked-read data from haplotagging, TELLseq, stLFR, and
even regular non-linked WGS data. It uses all the magic of [Snakemake](https://snakemake.readthedocs.io/en/stable/)
under the hood to handle  the worklfow decision-making, but as a user, you just interact with it like a normal command-line 
program. Harpy employs both well known and niche programs to take raw linked-read sequences and process
them to become called SNP genotypes (or haplotypes) or large structural variants (inversions, deletions, duplications).
Most of the settings are pre-configured and the settings you can modify are done at the command line. Some parts of this documentation
will refer to haplotagging specifically as we either forgot to update parts of the documentation or require you (the user)
to do a data conversion for some parts of Harpy to work with non-haplotagging linked-read data. As always, feel free to drop
an [Issue](https://github.com/pdimens/harpy/issues/new/choose) or open a [Discussion](https://github.com/pdimens/harpy/discussions) on GitHub.

## Harpy Commands
Harpy is modular, meaning you can use different parts of it independent from each other. Need to only align reads?
Great! Only want to call variants? Awesome! All modules are called by `harpy <workflow>`. For example, use `harpy align` to align reads.

{.compact}
| Command                                                                | Description                                                          |
|:------------------------------------------------------------------------|:---------------------------------------------------------------------|
| [!badge corners="pill" text="align"](Workflows/Align/Align.md)          | Align sample sequences to a reference genome                         |
| [!badge corners="pill" text="assembly"](Workflows/assembly.md)          | Create a genome assembly from linked-reads                           |
| [!badge corners="pill" text="convert"](Workflows/convert.md)            | Convert data between linked-read types                               |
| [!badge corners="pill" text="deconvolve"](Workflows/deconvolve.md)      | Resolve barcode sharing in unrelated molecules                       |
| [!badge corners="pill" text="downsample"](Workflows/downsample.md)      | Downsample data by barcode                                           |
| [!badge corners="pill" text="demultiplex"](Workflows/demultiplex.md)    | Demultiplex haplotagged FASTQ files                                  |
| [!badge corners="pill" text="impute"](Workflows/impute.md)              | Impute genotypes using variants and sequences                        |
| [!badge corners="pill" text="metassembly"](Workflows/metassembly.md)    | Create a metagenome assembly from linked-reads                       |
| [!badge corners="pill" text="phase"](Workflows/phase.md)                | Phase SNPs into haplotypes                                           |
| [!badge corners="pill" text="validate"](Workflows/validate.md)          | Run various format checks for FASTQ and BAM files                    |
| [!badge corners="pill" text="qc"](Workflows/qc.md)                      | Remove adapters, deduplicate, and quality trim sequences             |
| [!badge corners="pill" text="simulate"](Workflows/Simulate/Simulate.md) | Simulate linked reads or genomic variants                            |
| [!badge corners="pill" text="snp"](Workflows/snp.md)                    | Call SNPs and small indels                                           |
| [!badge corners="pill" text="sv"](Workflows/SV/SV.md)                   | Call large structural variants (inversions, deletions, duplications) |

## Using Harpy
You can call `harpy` without any arguments (or with `--help`) to print the docstring to your terminal. You can likewise call any of the modules without arguments or with `--help` to see their usage  (e.g. `harpy align --help`).
``` harpy --help                                                      
Usage: harpy COMMAND [ARGS]...

An automated workflow for linked-read data to go from raw data to
genotypes (or phased haplotypes). Batteries included.
demultiplex >> qc >> align >> snp >> impute >> phase >> sv
                                                                 
Documentation: https://pdimens.github.io/harpy/

Data Processing:
  align        Align sequences to a reference genome
  assembly     Assemble linked reads into a genome
  demultiplex  Demultiplex haplotagged FASTQ files
  impute       Impute variant genotypes from alignments
  metassembly  Assemble linked reads into a metagenome
  phase        Phase SNPs into haplotypes
  qc           FASTQ adapter removal, quality filtering, etc.
  simulate     Simulate genomic variants
  snp          Call SNPs and small indels from alignments
  sv           Call inversions, deletions, and duplications from alignments

Other Commands:
  deconvolve  Resolve barcode sharing in unrelated molecules
  template    Create files and HPC configs for workflows

Troubleshoot:
  deps      Locally install workflow dependencies
  diagnose  Attempt to resolve workflow errors
  resume    Continue an incomplete Harpy workflow
  validate  File format checks for linked-read data
  view      View a workflow's components
```

## Typical Linked-Read Workflows
Depending on your project goals, you may want any combination of SNPs, structural
variants (inversions, deletions, duplications), or phased haplotypes. Below are diagrams
outlining general workflows for linked-read data, depending on your goals.

+++ Variant Calling

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
+++ Assembly

```mermaid
graph LR
    QC([QC, trim adapters, etc.]):::clean--->DC([barcode deconvolution]):::clean
    DC--->Assembly([assembly/metassembly]):::clean

    classDef clean fill:#f5f6f9,stroke:#b7c9ef,stroke-width:2px
```
+++