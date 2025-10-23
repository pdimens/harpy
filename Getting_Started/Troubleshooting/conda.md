---
label: Install / Conda
icon: package
---

# :icon-package: Install / Conda

## Installation issues
Conda is an awesome package manager, but _was_ slow and used a ton of memory
as dependencies increased. Recent (`23.10+`) versions of Conda [now use the `libmamba` solver](https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community),
the super-fast and super-lightweight solver from Mamba. If you're experiencing
a suspiciously slow Harpy installation, either update Conda to at least version `23.10` or use Mamba/Pixi. If `conda` is forbidden
on your system, [see below](#conda-is-forbidden).

### strict channel priorities
Sometimes conda channel priorities may give you grief, be it when trying to install Harpy, or Harpy trying to install a particular
workflow's dependencies. Although Conda recommends setting channel priorities to `strict`, we find that it sometimes causes the solver
(within or without Snakemake) to fail to find the necessary software. One typical way to circumvent that is to set the channel priority to 
something other than `strict`:
```bash
conda config --set channel_priority true
# or #
conda config --set channel_priority false
```

## Conda issues
Unless using `--container`, Harpy leans on Snakemake to install small conda environments necessary to accomplish the tasks in workflows.
The most common issue we tend to see is errors during this process.

### installation into base
It's tempting to install software into the `base` environment conda/mamba ships with, but that behavior is discouraged (by Conda itself, not just us).
Modificatons to `base` can mess with all sorts of things in other environments ([see this issue](https://github.com/pdimens/harpy/issues/248)), so we recommend installing Harpy (or anything) into literally any other environment. If you find conda things just aren't working
right and you know you have at some point modified the `base` environment, you can reset your `base` environment. It's not guaranteed, but 
doing that might just fix things:
```bash
conda install --rev 0 --name base
```

## Install dependencies manually
### conda is forbidden
Some HPC configurations prohibit the use of `conda`. If that's preventing you from installing Harpy, then you can try to 
[install it as a python package](/Getting_Started/install.md#pip) using pip. You will still need to manually install the
package dependencies provided below. Since `conda` is disabled, you will need to use the `--container` option for most of
the commands.

==- additional harpy dependencies
{.compact}
| program             | minimum version | notes                     |
|:--------------------|:---------------:|:--------------------------|
| `bcftools`          |      1.22       | usually already installed |
| `htslib`            |      1.22       | usually already installed |
| `pysam`             |      0.23       |                           |
| `rich-click`        |       1.8       |                           |
| `snakemake-minimal` |       9.0       | `snakemake` works too     |
| `samtools`          |      1.22       | usually already installed |
| `seqtk`             |                 |                           |
===

### nodes do not have internet access
A not-uncommon HPC setup is to have the login node connected to the internet, but the worker nodes have no
internet access. This obviously creates a problem if Snakemake needs to download software dependencies (conda or container). To
facilitate this kind of setup, Harpy has the [!badge corners="pill" text="deps"]() command, which runs a fake workflow to download and
install workflow dependencies. The [!badge corners="pill" text="deps conda"]() and [!badge corners="pill" text="deps container"]() commands 
are meant to be used on the internet-connected login node in circumstances like this.

#### deps conda
Create the conda environments required by Harpy's workflows (e.g. `phase`).
Only useful for specific HPC configurations where worker nodes
do not have internet access to let snakemake install conda packages itself. 

```bash usage
harpy deps conda workflows...
```

```bash example | install all possible harpy workflow dependencies
harpy deps conda all
```

You can provide any combination of workflow names to selectively install only the necessary environments, or `all`
to install (you guessed it!) all of them: 
- `all`
- `align`
- `assembly`
- `metassembly`
- `phase`
- `qc`
- `r`
- `simulations`
- `stitch`
- `variants`

#### deps container
Manually pull the Harpy dependency container from dockerhub and convert it
into an Apptainer `.sif` file. Only useful for specific HPC configurations where worker nodes
do not have internet access to let snakemake download the container itself. 
Run this command again without arguments to use it.


```bash usage
harpy deps container
```