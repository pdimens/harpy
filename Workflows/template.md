---
label: Template
icon: repo-clone
description: Generate extra files for analysis with Harpy
order: 7
---

Creating template/boilerplate files has been consolidated into `harpy template`. These commands will
generate a specific file and write to `stdout`.

## groupings
Creates a sample grouping file for variant calling

```bash usage
harpy template groupings INPUTS > output
```

```bash example | generate a grouping file from folder of fastq/bam files
harpy template groupings data/ > samples.groups
```
### arguments
{.compact}
| argument   | description                                                                                   |
| :--------- | :-------------------------------------------------------------------------------------------- |
| `INPUTS`   | [!badge variant="info" text="required"] Files or directories containing input FASTQ/BAM files |

This optional file is useful if you want SNP variant calling to happen on a
per-population level via [!badge corners="pill" text="harpy snp"](snp.md#populations) or on samples
pooled-as-populations via [!badge corners="pill" text="harpy sv"](SV/naibr.md#pooled-sample-variant-calling).
- takes the format of sample[!badge variant="ghost" text="tab"]group
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
---
## hpc-*
Create template configurations for HPC cluster job submission systems (e.g. SLURM, HTConder) that can
be provided to the `--hpc` option for workflows. You will likely also need to install the appropriate
snakemake plugin(s) to use this feature. For example, to use the SLURM plugin, you will need to install
[`snakemake-executor-plugin-slurm`](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)
and likely [`snakemake-storage-plugin-fs`](https://snakemake.github.io/snakemake-plugin-catalog/plugins/storage/fs.html)
(via pip, conda, etc.).

```bash usage
harpy template hpc-system > out.yaml
```

```bash example | create SLURM submission template
harpy template hpc-slurm > slurm.yaml
```
---
## impute
Create a template parameter file for the [!badge corners="pill" text="impute"](/Workflows/impute.md) module. 
The file is formatted correctly and serves as a starting point for using parameters that make sense for your study.
Typically, one runs STITCH multiple times, exploring how results vary with
different model parameters. The solution Harpy uses for this is to have the user
provide a tab-delimited dataframe file where the columns are the 6 STITCH model 
parameters and the rows are the values for those parameters. To make formatting
easier, a template file is generated for you, just replace the values and add/remove
rows as necessary. See the section for the [!badge corners="pill" text="impute"](/Workflows/impute.md)
module for details on these parameters.

```bash usage
harpy template impute > output
```

```bash example | create imputation parameter template
harpy template impute > params.stitch
```

```text resulting params.stitch file
name model	usebx	bxlimit	k	s	ngen
k10_ng50 diploid	TRUE	50000	3	2	10
k1_ng30 diploid	TRUE	50000	3	1	5
high_ngen   diploid TRUE    50000   15  1   100
```

