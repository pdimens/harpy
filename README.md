[![logo](https://github.com/pdimens/harpy/blob/docs/static/logo_trans.png?raw=true)](https://pdimens.github.io/harpy)

[![GitHub Release](https://img.shields.io/github/v/release/pdimens/harpy?style=for-the-badge&logo=anaconda&logoColor=ffffff)](https://github.com/pdimens/harpy/releases)
[![documentation badge](https://img.shields.io/badge/read%20the-docs-fbab3a?style=for-the-badge&logo=quicklook&logoColor=ffffff)](https://pdimens.github.io/harpy)
[![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/pdimens/harpy/tests.yml?style=for-the-badge&logo=cachet&logoColor=ffffff)](https://www.youtube.com/watch?v=F1qdBPlK9M4)
[![haplotagging badge](https://custom-icon-badges.demolab.com/badge/-Haplotagging-8879b9?style=for-the-badge&logo=grapheneos&logoColor=ffffff)](https://www.fml.tuebingen.mpg.de/9418/haplotagging)

[Linked-read](https://doi.org/10.1073/pnas.2015005118) data processing pipeline. Getting you from raw linked reads to assemblies, genotypes, or phased haplotypes. Batteries included 🔋

<p align="center">
✨ Now works with TELLseq, stLFR, and non-linked read data! ✨
</p>

## 📥 Install
Detailed installation instructions are described in [the documentation](https://pdimens.github.io/harpy/install/). 
### 🐍 Conda
```bash
conda create -n harpy -c bioconda -c conda-forge harpy
conda activate harpy
```

### 🌟 Pixi
```bash
pixi global install -c conda-forge -c bioconda harpy
# or locally #
pixi init -c conda-forge -c bioconda projectname && cd projectname && pixi add harpy
```

### 📦 Containers
#### 🐳 Docker
```bash
docker pull quay.io/biocontainers/harpy
```

#### 🅰️ Apptainer
Find the most recent tag [here](https://quay.io/repository/biocontainers/harpy?tab=tags) and replace `$TAG` with it or use the `TAG=$(curl ...)` part below to pull the latest version name using the repository's API.
```bash
TAG=$(curl -s "https://quay.io/api/v1/repository/biocontainers/harpy/tag/" | cut -d'"' -f6)
apptainer pull docker://quay.io/biocontainers/harpy:$TAG
```

## ⚡ Usage
Just call `harpy` or `harpy --help` on the command line to get started! If installed via container, then call the containerized-Harpy however you are used to using containers on your system.

```bash
harpy module options... args...
```

## 🌈 Getting Started
No data? No problem! Use [HACk](https://davidebolo1993.github.io/visordoc/usecases/usecases.html#visor-hack) to simulate genomic variants from an existing genome and use [Mimick]((https://pdimens.github.io/mimick/#/)) to create linked-read data from an existing genome! You can see what haplotagging (or other linked read) data and Harpy are like without investing a single cent! A real-world walkthrough of how we did this for a benchmarking experiment can be found [here](https://pdimens.github.io/LRInversionSimulations/).

## Citation
>Dimens PV, Franckowiak RP, Iqbal A, Grenier JK, Munn PR, Therkildsen NO. Harpy: a pipeline for processing haplotagging linked-read data. Bioinform Adv. 2025 Jun 5;5(1):vbaf133. doi: [10.1093/bioadv/vbaf133](https://pubmed.ncbi.nlm.nih.gov/40575478/). PMID: 40575478; PMCID: PMC12198493.