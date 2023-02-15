![logo](https://github.com/pdimens/HARPY/blob/docs/static/logo.png?raw=true)

[![documentation badge](https://img.shields.io/badge/read%20the-documentation-c0bad4?style=for-the-badge&logo=Read%20The%20Docs)](https://pdimens.github.io/HARPY/#/) 

Experimental Haplotagging Data Processing Pipeline

#### Dependencies
Until this pipeline gets completed and hosted on Bioconda, it will be available by cloning/downloading this repository. The dependencies can be installed into a conda environment using the provided `harpyenv.yaml`:
```bash
conda env create --name harpy --file misc/harpyenv.yaml
```

The version of [EMA](https://github.com/arshajii/ema) bundled in this repository (`ema-h`) is a [fork](https://github.com/EdHarry/ema/tree/haplotag) of the orignal EMA modified to work with Generation 1 haplotag beadtags (AxxCxxBxxDxx).

#### Usage
```
./harpy --help
                                                           
 Usage: harpy [OPTIONS] COMMAND [ARGS]...                     
                                                              
 HARPY Haplotagging data processing pipeline.                 
 The pipeline trims reads, map sequences, calls variants,     
 imputes genotypes, and phases haplotypes.                   
                                                              
 trim 🡒 align 🡒 variants 🡒 impute 🡒 phase           
                                                              
 Documentation: https://pdimens.github.io/HARPY/#/                      
                                                              
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
