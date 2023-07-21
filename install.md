---
label: Install
icon: desktop-download
order: 100
---

# :icon-desktop-download: Install HARPY
=== :icon-checklist: You will need
A working installation of [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). 
We recommend you use [⚡mamba⚡](https://mamba.readthedocs.io/en/latest/installation.html). It's **so** much faster and uses a **lot** less memory. If using mamba, replace `conda` with `mamba` in the instructions below.
===

Harpy is now hosted on [Bioconda](https://anaconda.org/bioconda/harpy)! That means to install it, you just need to have `conda` (or `mamba`) on your Linux-based 
system and install it with a simple command. You can install Harpy into an existing environment or create a new one for it (recommended).

=== install into a new environment (**recommended**)
The code snippet below creates a new environment called `harpy` (the `-n harpy` part) and installs harpy into it from the bioconda channel (`-c bioconda` part). You can name this
environment anything (e.g. `haplotagging`, `jeanclaudevandamme`, etc.).
```bash install harpy
conda create -n harpy -c bioconda -c conda-forge harpy
```
Once conda finishes the install, you can activate the environment (with whatever `-n <env_name>` you gave it) and `harpy` will be callable
from the command line.
```bash activate harpy environment
conda activate <env_name>
```

==- install into an existing evironment
If you want to install harpy into an existing environment, then with an environment already activated (via `conda activate <env_name>`) simply use the `conda install` command and harpy
will be callable from the command line.

```bash install harpy
conda install -c bioconda harpy
```
==- local install for development
⚠️ **Not intended for regular users** ⚠️

If intent on installing harpy for development, you can do so by cloning the harpy
repository, installing the preconfigured conda environment, and running the `misc/buildlocal.sh`
script to move all the necessary files to the `/bin/` path within your active conda environment.

**Clone the repository**

```bash clone the repository
git clone https://github.com/pdimens/harpy.git
```
```bash install the dependencies with conda
conda env create --name harpy --file misc/harpyenv.yaml
```
This will create a conda environment named `harpy` with all the bits necessary to successfully run Harpy. You can change the name of this environment by specifying
`--name something`. 

**Activate the environment**

The environment with all the preinstalled dependencies can be activated with:
```bash activate the conda environment
# assuming the environment name is harpy from the step above
conda activate harpy
```
**Install the Harpy files**

Call the `misc/buildlocal.sh` bash script to install the Harpy-specific files into the conda environment:

```bash install the repo into the environment
bash misc/buildlocal.sh
```


===