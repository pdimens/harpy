---
label: Other
order: 1
icon: file-diff
description: Generate extra files for analysis with Harpy
---

# :icon-file-diff: Other Harpy modules
Some parts of Harpy (variant calling, imputation) want or need extra files. You can create various files necessary for different modules using these extra modules:

## :icon-terminal: Other modules
{.compact}
| module         | description                                                                      |
|:---------------|:---------------------------------------------------------------------------------|
| `popgroup`     | Create generic sample-group file using existing sample file names (fq.gz or bam) |
| `stitchparams` | Create template STITCH parameter file                                            |

### popgroup
Creates a sample grouping file for variant calling

```bash usage
harpy popgroup -o OUTPUTFILE INPUTS
```

```bash usage example
harpy popgroup -o samples.groups data/
```
#### arguments
{.compact}
| argument              | short name | type            | default | required | description                                                          |
|:----------------------|:----------:|:----------------|:-------:|:--------:|:---------------------------------------------------------------------|
| `INPUTS`           |            | file/directory paths  |         | **yes**  | Files or directories containing input FASTQ/BAM files             |
| `--output`               |    `-o`    | file path       |         | **yes**  | name of the output file                                             |

This optional file is useful if you want SNP variant calling to happen on a
per-population level via  [!badge corners="pill" text="harpy snp"](snp.md/#populations) or on samples
pooled-as-populations via [!badge corners="pill" text="harpy sv"](SV/naibr.md/#pooled-sample-variant-calling).
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
### stitchparams
Create a template parameter file for the [!badge corners="pill" text="impute"](/Modules/impute.md) module. The file is formatted correctly and serves
as a starting point for using parameters that make sense for your study.

```bash usage
harpy stitchparams -o OUTPUTFILE
```

```bash example
harpy stitchparams -o params.stitch
```

#### arguments
{.compact}
| argument              | short name | type            | default | required | description                                                          |
|:----------------------|:----------:|:----------------|:-------:|:--------:|:---------------------------------------------------------------------|
| `--output`               |    `-o`    | file path       |         | **yes**  | name of the output file                                             |

Typically, one runs STITCH multiple times, exploring how results vary with
different model parameters. The solution Harpy uses for this is to have the user
provide a tab-delimited dataframe file where the columns are the 6 STITCH model 
parameters and the rows are the values for those parameters. To make formatting
easier, a template file is generated for you, just replace the values and add/remove
rows as necessary. See the section for the [!badge corners="pill" text="impute"](/Modules/impute.md)
module for details on these parameters. The template file will look like:

``` params.stitch
model	usebx	bxlimit	k	s	ngen
diploid	TRUE	50000	3	2	10
diploid	TRUE	50000	3	1	5
```