---
label: Install
icon: desktop-download
order: 100
---

## :icon-desktop-download: Install Harpy
=== :icon-checklist: You will need
- [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [mamba](https://mamba.readthedocs.io/en/latest/installation.html)
  - we [strongly recommend](troubleshooting#installation-issue) you use conda `23.10` (or later) or mamba 
  - if using mamba, replace `conda` with `mamba` in the instructions below
===

Harpy is hosted on [Bioconda](https://anaconda.org/bioconda/harpy)! That means to install it, you just need to have `conda` (or `mamba`) on your Linux-based 
system and install it with a simple command. You can install Harpy into an existing environment or create a new one for it (recommended).

### install into a new environment 
##### ✨recommended✨
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

### install into an existing evironment
If you want to install harpy into an existing environment, then with an environment already activated (via `conda activate <env_name>`) simply use the `conda install` command and harpy
will be callable from the command line.

```bash install harpy
conda install -c bioconda -c conda-forge harpy
```
---

## :icon-move-to-top: Update Harpy
=== :icon-checklist: You will need
- [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [mamba](https://mamba.readthedocs.io/en/latest/installation.html)
- Harpy installed into an existing conda enviroment
===

If you want to update Harpy, the process is quite similar:
```bash update harpy
conda activate <env_name>
conda update -c conda-forge bioconda::harpy
```
