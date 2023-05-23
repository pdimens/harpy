[![logo](https://github.com/pdimens/harpy/blob/docs/static/logo.png?raw=true)](https://pdimens.github.io/harpy)

[![documentation badge](https://img.shields.io/badge/read%20the-documentation-fbab3a?style=for-the-badge&logo=Read%20The%20Docs)](https://pdimens.github.io/harpy) 

Experimental Haplotagging Data Processing Pipeline

## Install
Now hosted on [Bioconda](https://anaconda.org/bioconda/harpy)! Install harpy on Linux-based systems using [conda](https://mamba.readthedocs.io/en/latest/installation.html).

### into existing conda environment
If you wish to install harpy and its dependencies into an existing environment, activate that environment (`conda activate <env_name>`) and execute this `conda install` code:
```bash
conda install -c bioconda harpy

# or ⚡ mamba ⚡#
mamba install -c bioconda -c conda-forge harpy
```

### into new conda environment
To avoid dependency conflicts with an existing environment, it is best to create a new environment for a harpy installation. The code below creates a new conda environment called `harpy` (`-n harpy`) and installs harpy into it. You can name this environment whatever you like. 
```bash
conda create -n harpy -c bioconda harpy

# or ⚡ mamba ⚡#
mamba create -n harpy -c bioconda -c conda-forge harpy
```

### Activate the harpy environment
Once installed with one of the methods above, activate the conda environment you installed harpy into with
```bash
conda activate <env_name>
```
where `<env_name>` is the name of that environment. After doing so, the `harpy` executable should be callable from your path.

#### EMA note
The version of [EMA](https://github.com/arshajii/ema) bundled in this repository (`ema-h`) is a [fork](https://github.com/EdHarry/ema/tree/haplotag) of the orignal EMA modified to work with Generation 1 haplotag beadtags (AxxCxxBxxDxx). Work is underway to merge haplotag support and publish a new version of EMA to remove reliance on this precompiled fork. 

## Usage
```
harpy --help
                                                           
 Usage: harpy [OPTIONS] COMMAND [ARGS]...                     
                                                              
 HARPY Haplotagging data processing pipeline.                 
 The pipeline trims reads, map sequences, calls variants,     
 imputes genotypes, and phases haplotypes.                   
                                                              
 trim 🡒 align 🡒 variants 🡒 impute 🡒 phase           
                                                              
 Documentation: https://pdimens.github.io/harpy/                      
                                                              
╭─ Options ──────────────────────────────────────────────╮
│ --help      Show this message and exit.                │
╰────────────────────────────────────────────────────────╯
╭─ Commands ─────────────────────────────────────────────╮
│ align     Align sample sequences to a reference genome │
│ impute    Impute genotypes from genotype likelihoods   │
│ phase     Phase SNPs into haplotypes                   │
│ trim      Remove adapters and quality trim sequences   │
│ variants  Call variants from sample alignments         │
╰────────────────────────────────────────────────────────╯
```
