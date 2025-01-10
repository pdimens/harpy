[![logo](https://github.com/pdimens/harpy/blob/docs/static/logo_trans.png?raw=true)](https://pdimens.github.io/harpy)

[![GitHub Release](https://img.shields.io/github/v/release/pdimens/harpy?style=for-the-badge&logo=hackthebox&logoColor=ffffff)](https://github.com/pdimens/harpy/releases)
[![documentation badge](https://img.shields.io/badge/read%20the-documentation-fbab3a?style=for-the-badge&logo=searxng&logoColor=ffffff)](https://pdimens.github.io/harpy) 
[![Anaconda-Server Badge](https://img.shields.io/conda/dn/bioconda/harpy.svg?style=for-the-badge&logo=docusign&logoColor=ffffff)](https://anaconda.org/bioconda/harpy)
[![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/pdimens/harpy/tests.yml?style=for-the-badge&logo=cachet&logoColor=ffffff)](https://www.youtube.com/watch?v=F1qdBPlK9M4)

[Haplotagging](https://doi.org/10.1073/pnas.2015005118) Data Processing Pipeline. Getting you from raw reads to assemblies, genotypes, or phased haplotypes or your 💰 back.


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
install Harpy to be accessible in your PATH-- just make sure `~/.pixi/bin` is in your PATH:
```
# ~/.zshrc or ~/.bashrc (or equivalent) 
export PATH=~/.pixi/bin:$PATH
```
```bash
pixi global install -c conda-forge -c bioconda harpy
```

<details>
  <summary>⬇️ install as local environment </summary>

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

## ⚡ Usage
Just call `harpy` or `harpy --help` on the command line to get started!
```bash
harpy
```

## 🌈 Getting Started
No data? No problem! Harpy lets you [simulate genomic variants](https://pdimens.github.io/harpy/workflows/simulate/simulate-variants/)
from an existing genome and can also [create haplotag data](https://pdimens.github.io/harpy/workflows/simulate/simulate-linkedreads/)
from an existing genome! You can see what haplotag data (and Harpy) are like without paying a cent! A simple tutorial on simulating
both of these can be found [here](https://pdimens.github.io/harpy/blog/simulate_diploid/).
