---
label: Install
icon: desktop-download
order: 100
---

# :icon-desktop-download: Install Harpy
=== :icon-checklist: You will need one of either
{.compact}
| method                                                                                                                                                | considerations                                                                                                                                       |
|:------------------------------------------------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------|
| [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) or [mamba](https://mamba.readthedocs.io/en/latest/installation.html) | If using conda, we strongly recommend using version `23.10` or later. If using mamba, replace `conda` with `mamba` in the instructions below.        |
| [pixi](https://prefix.dev/blog/pixi_a_fast_conda_alternative)                                                                                         | Relatively new kid on the block-- think of it as the next-gen of conda. It's a little different, but it's _lightning fast_ and very sensible to use. |
| [docker](https://docs.docker.com/engine/install/)                                                                                                     | Clunkier to use, but best system compatibility and might be what your sysadmins prefer                                                               |

===

Harpy is hosted on [Bioconda](https://anaconda.org/bioconda/harpy), which means you just need to have either  `conda` or `pixi` (or `docker`) on your Unix-like 
system to install it.

+++ 🐍 Conda
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
Activate the conda environment with Harpy and call `harpy` or `harpy --help` on the command line to get started
```bash activate the environment
conda activate harpy_env

# or, if installed locally
conda activate path/to/harpy_env
```

```bash call harpy
harpy
```

+++ 🌟 Pixi

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
```bash
# global install
pixi global update harpy

# local install
# project dir has the pixi.toml file
cd path/to/projectdir
pixi update harpy
```
===

### ⚡ Pixi Usage
If installed globally, just call `harpy` or `harpy --help` on the command line, otherwise navigate to the directory with the pixi environment and call harpy from within it:
```bash call harpy if installed globally
# installed globally
harpy
```
```bash call harpy if installed locally
# navigate to project directory
cd path/to/workdir

# activate the pixi environment and run harpy
pixi shell
harpy

# or prepend `pixi run` to it without activating the environment
pixi run harpy
```

+++ 📦 Containers
If you didn't know, packages on Bioconda are automatically built as containers too!
### 🐳 Docker
```bash
docker pull quay.io/biocontainers/harpy
```

### 🅰️ Apptainer
Find the most recent tag [here](https://quay.io/repository/biocontainers/harpy?tab=tags) and replace `$TAG` with it or use the `TAG=$(curl ...)` part below to pull the latest version name using the repository's API.
```bash
TAG=$(curl -s "https://quay.io/api/v1/repository/biocontainers/harpy/tag/" | cut -d'"' -f6)
apptainer pull docker://quay.io/biocontainers/harpy:$TAG
```

Then proceed to use containerized-Harpy however you are used to using containers on your system.
+++