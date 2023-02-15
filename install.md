---
label: Install
icon: desktop-download
order: 100
---

# Install HARPY
||| :icon-checklist: You will need
A working installation of [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). We recommend you use [mamba](https://mamba.readthedocs.io/en/latest/installation.html). It's **so** much faster and uses a **lot** less memory.
|||

Until this pipeline gets completed and hosted on Bioconda, it will be available by cloning/downloading [the repository](https://github.com/pdimens/harpy). The dependencies can be installed into a conda environment using the provided `harpyenv.yaml`:
```bash
conda env create --name harpy --file misc/harpyenv.yaml
```

The version of [EMA](https://github.com/arshajii/ema) bundled in this repository (`ema-h`) is a [fork](https://github.com/EdHarry/ema/tree/haplotag) of the orignal EMA modified to work with Generation 1 haplotag beadtags (AxxCxxBxxDxx).