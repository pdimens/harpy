[![logo](https://github.com/pdimens/harpy/blob/docs/static/logo_trans.png?raw=true)](https://pdimens.github.io/harpy)

[![GitHub Release](https://img.shields.io/github/v/release/pdimens/harpy?style=for-the-badge&logo=hackthebox&logoColor=ffffff)](https://github.com/pdimens/harpy/releases)
[![documentation badge](https://img.shields.io/badge/read%20the-documentation-fbab3a?style=for-the-badge&logo=searxng&logoColor=ffffff)](https://pdimens.github.io/harpy) 
[![Anaconda-Server Badge](https://img.shields.io/conda/dn/bioconda/harpy.svg?style=for-the-badge&logo=docusign&logoColor=ffffff)](https://anaconda.org/bioconda/harpy)
[![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/pdimens/harpy/tests.yml?style=for-the-badge&logo=cachet&logoColor=ffffff)](https://www.youtube.com/watch?v=F1qdBPlK9M4)

[Haplotag](https://doi.org/10.1073/pnas.2015005118) Data Processing Pipeline. Getting you from raw reads to genotypes/phased haplotypes or your money back.


## 📥 Install 
To avoid dependency conflicts with an existing environment, it is best to create a new environment for a harpy installation. The code below creates a new conda/mamba environment called `harpy` (via `-n harpy`) and installs harpy into it. You can name this environment whatever you like using the `-n somename` argument. 
```bash
mamba create -n harpy -c bioconda -c conda-forge harpy
```

<details>
  <summary>⚪️ install into an existing conda environment ⚪️</summary>

  ---
  
If you wish to install harpy and its dependencies into an existing environment, activate that environment (`conda activate env_name`) and execute this installation code:
```bash
mamba install -c conda-forge bioconda::harpy
```
Or provide `-n envname` to install it into an existing environment named `envname`
```bash
mamba install -n envname -c conda-forge bioconda::harpy
```

---

</details>

## Update
```bash
mamba update -c conda-forge bioconda::harpy
```

## 🌟 Activate the harpy environment
Once conda/mamba finishes, activate the conda/mamba environment you installed harpy into with
```bash
conda activate env_name
```
where `env_name` is the name of that environment. After doing so, the `harpy` executable should be callable from your path.


## ⚡ Usage
Just call `harpy` or `harpy --help` on the command line to get started!
```bash
harpy
```

## 🌈 Get Started
No data? No problem! Harpy lets you [simulate genomic variants](https://pdimens.github.io/harpy/modules/simulate/simulate-variants/)
from an existing genome and can also [create haplotag data](https://pdimens.github.io/harpy/modules/simulate/simulate-linkedreads/)
from an existing genome! You can see what haplotag data (and Harpy) are like without paying a single cent!