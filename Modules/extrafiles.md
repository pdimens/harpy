---
order: 7
icon: file-diff
description: Generate extra files for analysis with Harpy
---

# :icon-file-diff: Generate Extra Files
Some parts of Harpy (variant calling, imputation) want or need extra files. You can create various files necessary for different modules using the `harpy extra` module:
```bash
harpy extra OPTIONS... 
```

The arguments represent different sub-commands and can be run in any order or combination to generate the files you need.

## :icon-terminal: Running Options
| argument          | short name | type           | default | required | description                                                                      |
|:------------------|:----------:|:---------------|:-------:|:--------:|:---------------------------------------------------------------------------------|
| `--popgroup`      |    `-p`    | folder path    |         |    no    | Create generic sample-group file using existing sample file names (fq.gz or bam) |
| `--stitch-params` |    `-s`    | file path      |         |    no    | Create template STITCH parameter file                                            |
| `--hpc`           |    `-h`    | string [slurm] |  slurm  |    no    | Create HPC scheduling profile for cluster submission                             |
| `--help`          |            |                |         |          | Show the module docstring                                                        |


### popgroup
||| `--popgroup`
**Sample grouping file for variant calling**

This file is entirely optional and useful if you want variant calling to happen on a per-population level using mpileup via `harpy variants -p`.
- takes the format of sample\<tab\>group
- all the samples will be assigned to group `1` since file names don't always provide grouping information, so make sure to edit the second column to reflect your data correctly.
- the file will look like:
```less popgroups.txt
sample1 1
sample2 1
sample3 2
sample4 1
sample5 3
```
|||

### stitch-params
||| `--stitch-params`
**STITCH parameter file**

Typically, one runs STITCH multiple times, exploring how results vary with
different model parameters. The solution Harpy uses for this is to have the user
provide a tab-delimited dataframe file where the columns are the 5 STITCH model 
parameters and the rows are the values for those parameters. To make formatting
easier, a template file is generated for you, just replace the values and add/remove
rows as necessary. See the [Imputation section](imputation.md) for details on these parameters.
|||

### hpc
||| `--hpc`
**HPC cluster profile**

For snakemake to work in harmony with an HPC scheduler, a "profile" needs to
be provided that tells Snakemake how it needs to interact with the HPC scheduler
to submit your jobs to the cluster. Using `harpy extra --hpc <hpc-type>` will create
the necessary folder and profile yaml file for you to use. To use the profile, call
the intended Harpy module with an extra ``--snakemake` argument:
```bash
# use the slurm profile
harpy module --option1 <value1> --option2 <value2> --snakemake "--profile slurm/"
```

#### Profiles
These are ongoing as I figure out how it works. Currently working on SLURM support, as that is
the system we have available.
- slurm: [by jdblischak](https://github.com/jdblischak/smk-simple-slurm)
|||