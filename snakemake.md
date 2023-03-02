---
label: Snakemake etc.
icon: gear
order: 2
---

# Adding additional Snakamake parameters
Harpy relies on Snakemake under the hood to handle file and job dependencies. Most of these details have been abstracted away from the end-user, but every module of Harpy (except `extra`) has an optional flag `-s` (`--snakemake`) that you can use to augment the Snakemake workflow if necessary. Whenever you use this flag, your argument must be enclosed in quotation marks, for example:
```bash
harpy trim -d rawseq -s "--dry-run"
```
This means you can add several Snakemake arguments at once, as long as the entire thing is enclosed in quotes:
```bash
harpy trim -d rawseq -s "--dry-run --debug --shadow-prefix /scratch"
```

## Reserved arguments
Harpy calls Snakemake using specific arguments, meaning you cannot append these again to the internal command line call. Well, you can, but Snakemake will error and exit. [Everything else](https://snakemake.readthedocs.io/en/stable/executing/cli.html#all-options) is allowed. The reserved (**forbidden**) arguments are:
- `--directory`
- `--cores`
- `--snakefile`
- `--config`

## Use cases
You likely wont need to invoke `--snakemake` very often, if ever. That being said, here are what might be the most common use cases for this parameter.

### Dry run
#### `--dry-run`
This is a directive in which Snakemake will build the DAG and "pretend" to run the Harpy workflow. Useful for knowing what you're getting yourself into ahead of time. It's also useful for debugging during development.

### Rerun an incomplete workflow
#### `--rerun-incomplete`
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