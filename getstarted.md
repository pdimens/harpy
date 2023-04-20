---
label: Before you start
icon: check-circle
order: 97
---
# :icon-check-circle: Before you start
## Required files and formats
Before you start using Harpy for your haplotagging data processing, you will need to make sure you have a few things on hand for things to go smoothly.
At the minimum, you will need:

==- 1. Haplotagging sequences, in b/gzipped FASTQ format
- file names must not have any wierd special characters or dot (`.`) separators in the sample name
    - legal: `sample_01_pop1.F.fq.gz`
    - legal: `sample-01.F.fq.gz`
    - **illegal**: `sample.01.pop1.F.fq.gz`
- the haplotagging sequences **must** have the barcode in the read headers. 
- the barcode must be in the format `AXXCXXBXXDXX`, where `XX` is a number between `00` and `96`
    - `00` indicates a missing/invalid barcode segment
- the barcode must be preceded by a `BX:Z:` tag in the read header
``` example header
@A00470:481:HNYFWDRX2:1:2101:16062:1031 BX:Z:A62C38B38D69 1:N:0:TATCAGTA+TTACTACT
```
==- 2. A reference genome, in FASTA format
A plain haploid genome assembly in uncompressed FASTA format.
===

## Adding additional Snakamake parameters
Harpy relies on Snakemake under the hood to handle file and job dependencies. Most of these details have been abstracted away from the end-user, but every module of Harpy (except `extra`) has an optional flag `-s` (`--snakemake`) that you can use to augment the Snakemake workflow if necessary. Whenever you use this flag, your argument must be enclosed in quotation marks, for example:
```bash
harpy trim -d rawseq -s "--dry-run"
```
This means you can add several Snakemake arguments at once, as long as the entire thing is enclosed in quotes:
```bash
harpy trim -d rawseq -s "--dry-run --debug --shadow-prefix /scratch"
```

!!!danger reserved/forbidden arguments
Harpy calls Snakemake with a given set of arguments, meaning you cannot append these again to the internal command line call. Well, you can, but Snakemake will error and exit. [Everything else](https://snakemake.readthedocs.io/en/stable/executing/cli.html#all-options) is allowed. The reserved (**forbidden**) arguments are:
- `--directory`
- `--cores`
- `--snakefile`
- `--config`
!!!

### Use cases
You likely wont need to invoke `--snakemake` very often, if ever. However, here are common use cases for this parameter.

==- Dry run
##### `--dry-run`
This is a directive in which Snakemake will build the DAG and "pretend" to run the Harpy workflow. Useful for knowing what you're getting yourself into ahead of time. It's also useful for debugging during development. Here is an example of dry-running variant calling:
```bash
harpy variants -g genome.fasta  -d ReadMapping/ema -s "--dry-run"
```
==- Rerun an incomplete workflow
##### `--rerun-incomplete`
There will be plenty of reasons that Harpy/Snakemake might end prematurely, like corrupt files, system errors, insufficient resources, etc.
When this happens, Snakemake has a save-state in the `.snakemake` folder where it knows the last run was incomplete, and when you try to run
Harpy again, Snakemake will complain. Something like this:
```
$ harpy variants --leviathan -g genome.fasta  -d ReadMapping/ema --threads 8 -p samples.groups
Building DAG of jobs...
IncompleteFilesException:
The files below seem to be incomplete. If you are sure that certain files are n
ot incomplete, mark them as complete with

    snakemake --cleanup-metadata <filenames>

To re-generate the files rerun your command with the --rerun-incomplete flag.
Incomplete files:
Variants/leviathan-pop/lrezIndexed/2.bci
```
So, the easiest workaround would be to regenerate the incomplete files and use `-s "--rerun-incomplete"`, like so:
```bash
harpy variants --leviathan -g genome.fasta  -d ReadMapping/ema --threads 8 -p samples.groups -s "--rerun-incomplete"
```
==- Specific file target
Sometimes you want to generate a specific intermediate file (or files) rather than running the entire module to completion. For example,
you want the beadtag report Harpy makes from the output of `EMA count`. To do this, just list the file/files (relative
to your working directory) without any flags. Example for the beadtag report:
```bash
harpy align -g genome.fasta -d Trimming -t 4 -s "Alignemnts/ema/stats/reads.bxstats.html"
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
harpy variants --leviathan -g genome.fasta  -d ReadMapping/ema --threads 8 -p samples.groups -s "--shadow-prefix /SCRATCH/username/"
```
==- Unlocking the directory
Sometimes Snakemake might scold/warn you about something you didn't realize you did. One
common case is when you prematurely terminate Harpy with `ctrl + c` or by terminating 
the process by other means. You might try to rerun Harpy afterwards and be met with a 
`LockException` like this:
```
Building DAG of jobs...
LockException:
Error: Directory cannot be locked. Please make sure that no other Snakemake pro
cess is trying to create the same files in the following directory:
/local/user/projectdir
If you are sure that no other instances of snakemake are running on this direct
ory, the remaining lock was likely caused by a kill signal or a power loss. It 
can be removed with the --unlock argument.
```
Like the error suggests, this can be overcome using the `--unlock` argument, which
would be provided to Harpy as `-s "--unlock"`. You may also remove
`.snakemake/` your working directory. Choose your own adventure, although unlocking
is likely safer. Unlocking the directory does not run Harpy, it just exists after unlocking (a snakemake nuance),
so you will need to run Harpy again  without this parameter.
===
