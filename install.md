---
label: Install
icon: desktop-download
order: 100
---

# :icon-desktop-download: Install HARPY
=== :icon-checklist: You will need
A working installation of [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html). 
We [strongly recommend](issues#problem-installing-with-conda) you use [⚡mamba⚡](https://mamba.readthedocs.io/en/latest/installation.html). It's **so** much faster and uses a **lot** less memory. If using conda, replace `mamba` with `conda` in the instructions below.
===

Harpy is now hosted on [Bioconda](https://anaconda.org/bioconda/harpy)! That means to install it, you just need to have `mamba` (or `conda`) on your Linux-based 
system and install it with a simple command. You can install Harpy into an existing environment or create a new one for it (recommended).

=== install into a new environment (**recommended**)
The code snippet below creates a new environment called `harpy` (the `-n harpy` part) and installs harpy into it from the bioconda channel (`-c bioconda` part). You can name this
environment anything (e.g. `haplotagging`, `jeanclaudevandamme`, etc.).
```bash install harpy
mamba create -n harpy -c bioconda -c conda-forge harpy
```
Once conda finishes the install, you can activate the environment (with whatever `-n <env_name>` you gave it) and `harpy` will be callable
from the command line.
```bash activate harpy environment
mamba activate <env_name>
```

==- install into an existing evironment
If you want to install harpy into an existing environment, then with an environment already activated (via `conda activate <env_name>`) simply use the `conda install` command and harpy
will be callable from the command line.

```bash install harpy
mamba install -c bioconda harpy
```
===