---
label: Install
icon: desktop-download
order: 100
---

# :icon-desktop-download: Install Harpy
=== :icon-checklist: You will need one of either
- [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [mamba](https://mamba.readthedocs.io/en/latest/installation.html)
  - if using conda, we strongly recommend using version `23.10` or later 
  - if using mamba, replace `conda` with `mamba` in the instructions below
- [pixi](https://prefix.dev/blog/pixi_a_fast_conda_alternative)
  - relatively new kid on the block-- think of it as the next-gen of conda
  - it's a little different, but it's _lightning fast_ and very sensible to use
===

Harpy is hosted on [Bioconda](https://anaconda.org/bioconda/harpy), which means you just need to have either  `conda` or `pixi` on your Linux-based 
system to install it. This page details both the conda and pixi installation approaches.

## 🐍 Conda
It's best to create a new environment for a harpy installation. The code below creates a new conda/mamba environment called `harpy` (via `-n harpy`) and installs harpy into it. You can name this environment whatever you like using the `-n somename` argument. 
```bash
conda create -n harpy -c bioconda -c conda-forge harpy
```

Once conda/mamba finishes, activate the harpy conda/mamba environment with:
```bash
conda activate env_name
```
where `env_name` is the name of that environment. After doing so, the `harpy` executable should be callable from your path.

==- ⬇️ install as local conda environment
Alternatively, you can create the environment locally within a specific project folder, just swap `-n harpy` for
`-p path/to/workdir/harpy`, which creates the environment in that specific folder (e.g. `potato_blight/harpy`).
```
conda create -p path/to/workdir/harpy -c bioconda -c conda-forge harpy
```
==- ⬇️ install into existing conda environment
 
If you wish to install harpy and its dependencies into an existing environment, activate that environment and execute this installation code:
```bash
conda activate env_name
conda install -c conda-forge bioconda::harpy
```
Or provide `-n env_name` to install it into an existing environment named `env_name`
```bash
conda install -n env_name -c conda-forge bioconda::harpy
```
==- ⬆️ updating harpy

Activate the environment and run `conda update`:

```bash
conda update -c conda-forge bioconda::harpy
```
===

### ⚡ Conda Usage
+++ Installed Globally

Activate the conda environment with Harpy and call `harpy` or `harpy --help` on the command line to get started
```bash activate the environment
conda activate harpy_env
```
```bash call harpy
harpy
```

+++ Installed Locally
Activate the conda environment with Harpy and call `harpy` or `harpy --help` on the command line to get started
```bash activate the environment
conda activate path/to/harpy_env
```

```bash call harpy
harpy
```
+++

----

## 🌟 Pixi
If you prefer [Pixi](https://pixi.sh/latest/) (it's pretty good, you should try it), you can
install Harpy to be accessible in your PATH, i.e. a "global" installation:

```bash
pixi global install -c conda-forge -c bioconda harpy
```
!!! add pixi to PATH
Make sure `~/.pixi/bin` is in your PATH:
```bash ~/.zshrc or ~/.bashrc (or equivalent) 
export PATH=~/.pixi/bin:$PATH
```
!!!

==- ⬇️ install as local environment

Likewise, you can do an installation into a local project directory:

```bash
pixi init -c conda-forge -c bioconda projectname && cd projectname
pixi add harpy
```

==- ⬆️ updating harpy

If installed via Pixi, you can update Harpy with `pixi update`:
+++ global install
```bash
pixi global update harpy
```
+++ local install
```bash
# project dir has the pixi.toml file
cd path/to/projectdir
pixi update harpy
```
+++
===

### ⚡ Pixi Usage
+++ Installed Globally
If installed globally, just call `harpy` or `harpy --help` on the command line:
```bash
harpy
```
+++ Installed Locally
```bash Navigate to project directory
cd path/to/workdir
```

Then:

```bash call harpy from within pixi environment
pixi shell

harpy --help
```
or
```bash call harpy with pixi prefix
pixi run harpy
```
+++