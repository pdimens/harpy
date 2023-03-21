---
label: Install
icon: desktop-download
order: 100
---

# Install HARPY
=== :icon-checklist: You will need
A working installation of [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). We recommend you use [mamba](https://mamba.readthedocs.io/en/latest/installation.html). It's **so** much faster and uses a **lot** less memory. If using mamba, replace `conda` with `mamba` in the instructions below.
===


## Clone the repository
Until this pipeline gets completed and hosted on Bioconda, it will be available by cloning/downloading [the repository](https://github.com/pdimens/harpy). 
```bash clone the repository
git clone https://github.com/pdimens/HARPY.git
```
## Install the dependencies
The dependencies can be installed into a conda environment using the provided `harpyenv.yaml`:
```bash install the dependencies with conda
conda env create --name harpy --file misc/harpyenv.yaml
```
This will create a conda environment named `harpy` with all the bits necessary to successfully run Harpy. You can change the name of this environment by specifying
`--name something`. 

## Activate the environment
The environment with all the preinstalled dependencies can be activated with:
```bash activate the conda environment
# assuming the environment name is harpy from the step above
conda activate harpy
```

!!!info EMA fork
The version of [EMA](https://github.com/arshajii/ema) bundled in this repository (`ema-h`) is a [fork](https://github.com/EdHarry/ema/tree/haplotag) of the orignal EMA modified to work with Generation 1 haplotag beadtags (AxxCxxBxxDxx). You will need to add the `ema-h` binary to your `PATH`. e.g. `cp ema/ema-h /usr/local/bin`.
!!!