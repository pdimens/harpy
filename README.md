[![logo](https://github.com/pdimens/harpy/blob/docs/static/logo_trans.png?raw=true)](https://pdimens.github.io/harpy)

![GitHub Release](https://img.shields.io/github/v/release/pdimens/harpy?style=for-the-badge)
[![documentation badge](https://img.shields.io/badge/read%20the-documentation-fbab3a?style=for-the-badge&logo=Read%20The%20Docs)](https://pdimens.github.io/harpy) 
[![Anaconda-Server Badge](https://img.shields.io/conda/dn/bioconda/harpy.svg?style=for-the-badge)](https://anaconda.org/bioconda/harpy)

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
mamba install -c bioconda -c conda-forge harpy
```
Or provide `-n envname` to install it into an existing environment named `envname`
```bash
mamba install -n envname -c bioconda -c conda-forge harpy
```

---

</details>

### 🌟 Activate the harpy environment
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
