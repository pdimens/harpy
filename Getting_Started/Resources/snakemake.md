---
label: Snakemake Things
icon: git-merge
---

# :icon-git-merge: Snakamake Things
## Workflow logs
Barring a few exceptions, most of Harpy's options are Snakemake workflows.
This means we are all at the mercy of how Snakemake operates, which includes
the `.snakemake/` folder in your project directory. That folder contains
all sorts of things necessary for Snakemake to do its magic. However, as a
convenience, Harpy workflows will move the Snakemake workflow log (all those things that print on screen when snakemake is running)
in a workflow's output directory. These logs are found in `OUTDIR/logs/snakemake`
and are named `workflow.X.DATE.log`, where `workflow` is the harpy workflow
(qc, sv_naibr, etc.), `X` is the attempt number (given by `X`, e.g. `4`), and
`DATE` is given as `MONTH_DAY_YEAR`. As an example, using the default
settings of `harpy qc`, you will find your workflow's log in
`QC/logs/snakemake/qc.1.5_13_2024.log.gz`. 

!!!ghost Design choices
Snakemake prints **a lot** of text during runtime. For some workflows, the resulting log files
could occupy >1 gigabyte of hard drive space 😱. For this reason, Harpy gzip compresses the log file
and deletes the original Snakemake logfile that you would expect to find in `.snakemake/logs/`.
!!!

## Adding Snakemake Parameters
Harpy relies on Snakemake under the hood to handle file and job dependencies.
Most of these details have been abstracted away from the end-user, but most
Harpy modules have an optional flag `--snakemake` that you can use to augment
the Snakemake workflow if necessary. [The full list of Snakemake command line
options can be found here](https://snakemake.readthedocs.io/en/stable/executing/cli.html).
Whenever you use this flag, your argument must be enclosed in quotation marks, for example:
```bash
harpy qc --snakemake "--dry-run" rawseq
```
This means you can add several Snakemake arguments at once, as long as the entire thing is enclosed in quotes:
```bash
harpy qc --snakemake "--dry-run --debug --shadow-prefix /scratch" rawseq
```

### Common use cases
You likely wont need to invoke `--snakemake` very often, if ever. However, 
here examples of some possible use cases for this parameter.

==- Specific file target
Sometimes you want to generate a specific intermediate file (or files) rather than running the entire module to completion. For example,
you want the beadtag report Harpy makes from the output of `EMA count`. To do this, just list the file/files (relative
to your working directory) without any flags. Example for the beadtag report:
```bash
harpy align bwa -g genome.fasta -t 4 --snakemake "Align/ema/reports/bxstats.html" QC/
```
This of course necessitates knowing the names of the files ahead of time. See the individual workflows for a breakdown of expected outputs. 

==- Set a shadow directory
##### `--shadow-prefix <dirname>`
If running Harpy on an HPC, your system administrator may enforce a policy that all data needs to be moved to a particular
network-attached storage for execution. On some systems this is called a `SCRATCH/` drive (or something similar) and files
are automatically deleted from that drive after the completion of an HPC-scheduled job. Rather than manually adding and removing
your files into that workspace storage to make sure you keep the output of your work while not running afoul of your HPC's policy,
you may use `--shadow-prefix <dirname>` where `<dirname>` is the path to the mandatory directory you need to work out of. By 
configuring this "shadow directory" setting, Snakemake will automatically move the files in/out of that directory for you:
```bash
harpy sv leviathan -g genome.fasta --threads 8 -p samples.groups --snakemake "--shadow-prefix /SCRATCH/username/" Align/bwa
```
===
