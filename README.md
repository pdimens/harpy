[![logo](https://github.com/pdimens/harpy/blob/docs/static/logo.png?raw=true)](https://pdimens.github.io/harpy)

[![documentation badge](https://img.shields.io/badge/read%20the-documentation-fbab3a?style=for-the-badge&logo=Read%20The%20Docs)](https://pdimens.github.io/harpy) 
[![Anaconda-Server Badge](https://img.shields.io/conda/dn/bioconda/harpy.svg?style=for-the-badge)](https://anaconda.org/bioconda/harpy)

[Haplotag](https://doi.org/10.1073/pnas.2015005118) Data Processing Pipeline. Getting you from raw reads to genotypes/phased haplotypes or your money back.

## Install
Now hosted on [Bioconda](https://anaconda.org/bioconda/harpy)! Install harpy on Linux-based systems using [conda](https://mamba.readthedocs.io/en/latest/installation.html) or [mamba](https://mamba.readthedocs.io/en/latest/micromamba-installation.html#umamba-install) (**recommended**).

### Install into new environment (recommended)
To avoid dependency conflicts with an existing environment, it is best to create a new environment for a harpy installation. The code below creates a new conda environment called `harpy` (via `-n harpy`) and installs harpy into it. You can name this environment whatever you like using the `-n somename` argument. 
```bash
conda create -n harpy -c bioconda harpy

# or ⚡ mamba ⚡
mamba create -n harpy -c bioconda -c conda-forge harpy
```

<details>
  <summary>install into an existing conda environment</summary>
 
### Install into existing environment
If you wish to install harpy and its dependencies into an existing environment, activate that environment (`conda activate <env_name>`) and execute this `conda install` code:
```bash
conda install -c bioconda harpy

# or ⚡ mamba ⚡
mamba install -c bioconda -c conda-forge harpy
```
</details>

### Activate the harpy environment
Once installed with one of the methods above, activate the conda environment you installed harpy into with
```bash
conda activate <env_name>
```
where `<env_name>` is the name of that environment. After doing so, the `harpy` executable should be callable from your path.

## Usage
Just call `harpy` or `harpy --help` on the command line to get started!

```                                                                 
 Usage: harpy COMMAND [ARGS]...                                  
                                                                 
                                                                 
                   Harpy haplotagging pipeline                   
 An automated workflow to demultiplex sequences, trim reads, 
 map sequences, call variants, impute genotypes, and phase 
 haplotypes of Haplotagging data. Batteries included.                          
                                                                 
 demultiplex >> trim >> align >> variants >> impute >> phase     
                                                                 
 Documentation: https://pdimens.github.io/harpy/                 
                                                                 
╭─ Options ─────────────────────────────────────────────────────╮
│ --version      Show the version and exit.                     │
│ --help     -h  Show this message and exit.                    │
╰───────────────────────────────────────────────────────────────╯
╭─ Commands ────────────────────────────────────────────────────╮
│ align        Align sample sequences to a reference genome     │
│ demultiplex  Demultiplex haplotagged FASTQ files              │
│ extra        Create various optional/necessary input files    │
│ impute       Impute genotypes using variants and sequences    │
│ phase        Phase SNPs into haplotypes                       │
│ trim         Remove adapters and quality trim sequences       │
│ variants     Call variants (SNP/SV) from samples              │
╰───────────────────────────────────────────────────────────────╯
```