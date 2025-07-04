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
### 🐍 Conda
It's best to create a new environment for a harpy installation. The code below creates a new conda/mamba environment called `harpy` (via `-n harpy`) and installs harpy into it. You can name this environment whatever you like using the `-n somename` argument. 
```bash
conda create -n harpy -c bioconda -c conda-forge harpy
```

Once conda/mamba finishes, activate the harpy conda/mamba environment with:
```bash
conda activate env_name
```
where `env_name` is the name of that environment. After doing so, the `harpy` executable should be callable from your path.

<details>
  <summary>⬇️ install as local conda environment </summary>

Alternatively, you can create the environment locally within a specific project folder, just swap `-n harpy` for
`-p path/to/workdir/harpy`, which creates the environment in that specific folder (e.g. `potato_blight/harpy`).
```
# for local project directory
conda create -p path/to/workdir/harpy -c bioconda -c conda-forge harpy
```

</details>

<details>
  <summary>⬇️ install into existing conda environment </summary>
 
If you wish to install harpy and its dependencies into an existing environment, activate that environment (`conda activate env_name`) and execute this installation code:
```bash
conda install -c conda-forge bioconda::harpy
```
Or provide `-n envname` to install it into an existing environment named `envname`
```bash
conda install -n envname -c conda-forge bioconda::harpy
```

</details>

<details>
  <summary>⬆️ updating harpy </summary>

If installed via conda, you can update Harpy by activating the environment
and running `conda update` like so:

```bash
conda update -c conda-forge bioconda::harpy
```

</details>

### 🌟 Pixi
If you prefer [Pixi](https://pixi.sh/latest/) (it's pretty good, you should try it), you can
install Harpy to be accessible in your PATH:

<details>
  <summary>🌟 how to install pixi </summary>

```bash
# install pixi
curl -fsSL https://pixi.sh/install.sh | bash

# add this to ~/.zshrc or ~/.bashrc (or equivalent) 
export PATH=~/.pixi/bin:$PATH
```

</details>


```bash
pixi global install -c conda-forge -c bioconda harpy
```

<details>
  <summary>⬇️ install harpy into local environment </summary>

Likewise, you can do an installation into a local project directory:

```bash
pixi init -c conda-forge -c bioconda projectname && cd projectname
pixi add harpy
```
After that finishes, you can activate the environment with:
```bash
pixi shell
```
Or run `harpy` by prefixing it with `pixi run`:
```bash
pixi run harpy
```
</details>

<details>
  <summary>⬆️ updating harpy </summary>

If installed via Pixi, you can update Harpy with `pixi update`:
```bash
# global install
pixi global update harpy

# local install
# project dir has the pixi.toml file
cd path/to/projectdir
pixi update harpy
```

</details>

### 📦 Docker
If you didn't know, packages on Bioconda are automatically built as containers too! So, if you're using docker, you can pull the Harpy container using:
```bash
docker pull quay.io/biocontainers/harpy
```
Then proceed to use containerized-Harpy however you are used to using containers on your system.

## ⚡ Usage
Just call `harpy` or `harpy --help` on the command line to get started!
```bash
harpy
```

## 🌈 Getting Started
No data? No problem! Harpy lets you [simulate genomic variants](https://pdimens.github.io/harpy/workflows/simulate/simulate-variants/)
from an existing genome and can also [create linked-read data](https://pdimens.github.io/harpy/workflows/simulate/simulate-linkedreads/)
from an existing genome! You can see what haplotagging (or other linked read) data and Harpy are like without paying a cent! A simple tutorial on simulating
both of these can be found [here](https://pdimens.github.io/harpy/blog/simulate_diploid/).
