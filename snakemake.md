---
label: Sneaky Snakemake
icon: terminal
order: 49
---

# :icon-terminal: Adding Snakamake parameters
Harpy relies on Snakemake under the hood to handle file and job dependencies.
Most of these details have been abstracted away from the end-user, but every
module of Harpy (except `extra`) has an optional flag `-s` (`--snakemake`) 
that you can use to augment the Snakemake workflow if necessary. Whenever you
use this flag, your argument must be enclosed in quotation marks, for example:
```bash
harpy trim -d rawseq -s "--dry-run"
```
This means you can add several Snakemake arguments at once, as long as the entire thing is enclosed in quotes:
```bash
harpy trim -d rawseq -s "--dry-run --debug --shadow-prefix /scratch"
```

!!!danger reserved/forbidden arguments
Harpy calls Snakemake with a given set of arguments, meaning you cannot append
these again to the internal command line call. Well, you can, but Snakemake will
error and exit. [Everything else](https://snakemake.readthedocs.io/en/stable/executing/cli.html#all-options) 
is allowed. The reserved (**forbidden**) arguments are:
- `--directory`
- `--cores`
- `--snakefile`
- `--config`
- `--rerun-incomplete`
- `--nolock`
!!!

### Common use cases
You likely wont need to invoke `--snakemake` very often, if ever. However, 
here are common use cases for this parameter.

==- Dry run
##### `--dry-run`
This is a directive in which Snakemake will build the DAG and "pretend" to
run the Harpy workflow. Useful for knowing what you're getting yourself into
ahead of time. It's also useful for debugging during development. Here is an 
example of dry-running variant calling:
```bash
harpy variants snp -g genome.fasta  -d Align/ema -s "--dry-run"
```
==- Specific file target
Sometimes you want to generate a specific intermediate file (or files) rather than running the entire module to completion. For example,
you want the beadtag report Harpy makes from the output of `EMA count`. To do this, just list the file/files (relative
to your working directory) without any flags. Example for the beadtag report:
```bash
harpy align -g genome.fasta -d Trim/ -t 4 -s "Alignemnts/ema/stats/reads.bxstats.html"
```
This of course necessitates knowing the names of the files ahead of time. See the individual modules for a breakdown of expected outputs. 

==- Set a shadow directory
##### `--shadow-prefix <dirname>`
If running Harpy on an HPC, your system administrator may enforce a policy that all data needs to be moved to a particular
network-attached storage for execution. On some systems this is called a `SCRATCH/` drive (or something similar) and files
are automatically deleted from that drive after the completion of an HPC-scheduled job. Rather than manually adding and removing
your files into that workspace storage to make sure you keep the output of your work while not running afoul of your HPC's policy,
you may use `--shadow-prefix <dirname>` where `<dirname>` is the path to the mandatory directory you need to work out of. By 
configuring this "shadow directory" setting, Snakemake will automatically move the files in/out of that directory for you:
```bash
harpy variants sv --method leviathan -g genome.fasta  -d Align/bwa --threads 8 -p samples.groups -s "--shadow-prefix /SCRATCH/username/"
```
===
