---
label: Other
order: 1
icon: file-diff
description: Generate extra files for analysis with Harpy
---

# :icon-file-diff: Other Harpy modules
Some parts of Harpy (variant calling, imputation) want or need extra files. You can create various files necessary for different modules using these extra modules:
The arguments represent different sub-commands and can be run in any order or combination to generate the files you need.

## :icon-terminal: Other modules
| module         | description                                                                      |
|:---------------|:---------------------------------------------------------------------------------|
| `popgroup`     | Create generic sample-group file using existing sample file names (fq.gz or bam) |
| `stitchparams` | Create template STITCH parameter file                                            |
| `hpc`          | Create HPC scheduling profile for cluster submission                             |

### popgroup
#### Sample grouping file for variant calling
##### arguments
- `-o`, `--output`: name of the output file
- `-d`, `--directory`: name of the directory of input files, either fastq or bam.

This file is entirely optional and useful if you want SNP variant calling to happen on a
per-population level via `harpy snp ... -p` or on samples pooled-as-populations via `harpy sv ... -p`.
- takes the format of sample\<tab\>group
- all the samples will be assigned to group `pop1` since file names don't always provide grouping information
    - so make sure to edit the second column to reflect your data correctly.
- the file will look like:
```less popgroups.txt
sample1 pop1
sample2 pop1
sample3 pop2
sample4 pop1
sample5 pop3
```

### stitchparams
#### STITCH parameter file
##### arguments
- `-o`, `--output`: name of the output file

Typically, one runs STITCH multiple times, exploring how results vary with
different model parameters. The solution Harpy uses for this is to have the user
provide a tab-delimited dataframe file where the columns are the 6 STITCH model 
parameters and the rows are the values for those parameters. To make formatting
easier, a template file is generated for you, just replace the values and add/remove
rows as necessary. See the [Imputation section](/Modules/impute.md) for details on these parameters.

### hpc
#### HPC cluster profile
!!!warning
HPC support is not yet natively integrated into Harpy. Until then, you can manually
use the [Snakemake HPC infrastructure](https://snakemake.readthedocs.io/en/stable/executing/cluster.html) with the `-s` flag.
!!!
##### arguments
- `-o`, `--output`: name of the output file
- `-s`, `--system`: name of the scheduling system
    - options: `slurm` (more to come)

For snakemake to work in harmony with an HPC scheduler, a "profile" needs to
be provided that tells Snakemake how it needs to interact with the HPC scheduler
to submit your jobs to the cluster. Using `harpy hpc -s <hpc-type>` will create
the necessary folder and profile yaml file for you to use. To use the profile, call
the intended Harpy module with an additional ``--snakemake` argument:
```bash use the slurm profile
harpy module --option1 <value1> --option2 <value2> --snakemake "--profile slurm.profile"
```