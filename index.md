# :icon-home: Home

![](static/logo_trans.png)

Harpy is a [linked-read](Getting_Started/inputformat.md) and non-linked WGS data processing pipeline for Linux-based systems.
It uses all the magic of [Snakemake](https://snakemake.readthedocs.io/en/stable/)
under the hood to handle the worklfow decision-making, but as a user, you just interact with it like a normal command-line 
program. Harpy employs both well known and niche programs to take raw linked-read sequences and process
them to become called SNP genotypes (or haplotypes) or large structural variants (inversions, deletions, duplications).
Feel free to open an [Issue](https://github.com/pdimens/harpy/issues/new/choose) or begin a [Discussion](https://github.com/pdimens/harpy/discussions) on GitHub.

[!card icon= ":heart:" title="Harpy is friendly" text="Drawing on the lessons of its predecessors and contemporaries, we strive to minimize the commonplace **struggle** of bioinformatics, inasmuch as we can." layout="compact"]("")
[!card icon= ":keyboard:" title="Harpy is hackable" text="Harpy's commands expose the most common and consequential arguments of the key software it will be running, but the workflows also 100% hackable." layout="compact"]("")
[!card icon= ":bento:" title="Harpy is modular" text="We believe in the 'pause and assess' approach between data processing steps, which is why there is such a robust report system." layout="compact"]("")
[!card icon= ":computer:" title="Harpy is not for analysis" text="Harpy leverages linked-read data to get you as far as genotypes or assemblies, without making assumptions about whether it will be used for popgen, biomed, etc." layout="compact"]("")

## Commands
Harpy is modular, meaning you can use different parts of it independent from each other. Need to only align reads?
Great! Only want to call variants? Awesome! All modules are called by `harpy <workflow>`. For example, use `harpy align` to align reads.
You can call `harpy` without any arguments (or with `--help`) to print the docstring to your terminal. You can likewise call any of the modules without arguments or with `--help` to see their usage, e.g.:
```bash
harpy align --help
```

## Utilities
An installation of Harpy also includes a series of [scripts/utilities](Getting_Started/Resources/utilities.md) called `harpy-utils` that are available along with the `harpy` package. These scripts are used within Harpy workflows, but you can also use them outside of Harpy workflows.
```bash
harpy-utils molecule-coverage
```

## Typical Workflows
Depending on your project goals, you may want any combination of SNPs, structural
variants (inversions, deletions, duplications), or phased haplotypes. Below are diagrams
outlining general workflows for linked-read data, depending on your goals.

+++ Variant Calling
>>> Preprocess
Sample demultiplexing and linked-read barcode demultiplexing

>>> QC
Remove adapters, low quality sequences, reads that are too short, poly-G tails, etc.

>>> Align
Align sequences to a reference genome

>>> SNP
Call Single Nucleotide Polymorphisms and small indels from alignments

>>> Impute (optional)
Use existing data to heuristically fill missing data

>>> Phase
Convert individual SNPs into multi-allele haplotypes reflecting alleles that were inherited together from each parent

>>> SV
Call structural variants (inversions, large deletions, and duplications) from alignments
>>>

+++ Assembly

>>> Preprocess
Sample demultiplexing and linked-read barcode demultiplexing

>>> QC
Remove adapters, low quality sequences, reads that are too short, poly-G tails, etc.

>>> Deconvolute
Correct linked-read barcodes for unrelated sequences that share the same barcode by chance ("clashing")

>>> Assembly
Assemble sequences into a genome or metagenome
>>>

+++
