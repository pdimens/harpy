---
label: Other
icon: ellipsis
description: Generate extra files for analysis with Harpy
order: 7
---

# :icon-ellipsis: Other Harpy modules
On this page you'll find Harpy functions that do other, ancillary things.

{.compact}
| module     | description                                              |
|:-----------|:---------------------------------------------------------|
| `diagnose` | Run the Snakemake debugger to identify hang-ups          |
| `resume`   | Continue a Harpy workflow from an existing output folder |
| `template` | Create template Harpy files                              |
| `view`     | View a workflow log, config, or snakefile                |


## diagnose
This will run dry-run Snakemake for a workflow with the `--dry-run --debug-dag` options
and print output to try to identify why Snakemake might be stalling for a workflow.
This is typically used during development or for troubleshooting edge-cases and you're likely/hopefully
not going to need to use this feature.

---

## template
Creating template/boilerplate files has been consolidated into `harpy template`. These commands will
generate a specific file and write to `stdout`.

### groupings
Creates a sample grouping file for variant calling

```bash usage
harpy template groupings INPUTS > output
```

```bash usage example
harpy template groupings data/ > samples.groups
```
#### arguments
{.compact}
| argument   | description                                                                                   |
| :--------- | :-------------------------------------------------------------------------------------------- |
| `INPUTS`   | [!badge variant="info" text="required"] Files or directories containing input FASTQ/BAM files |

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

### hpc-*
Create template configurations for HPC cluster job submission systems (e.g. SLURM, HTConder) that can
be provided to the `--hpc` option for workflows.

```bash usage
harpy template hpc-system > out.yaml
```

```bash example
harpy template hpc-slurm > slurm.yaml
```

### impute
Create a template parameter file for the [!badge corners="pill" text="impute"](/Workflows/impute.md) module. 
The file is formatted correctly and serves as a starting point for using parameters that make sense for your study.

```bash usage
harpy template impute > output
```

```bash example
harpy template impute > params.stitch
```

Typically, one runs STITCH multiple times, exploring how results vary with
different model parameters. The solution Harpy uses for this is to have the user
provide a tab-delimited dataframe file where the columns are the 6 STITCH model 
parameters and the rows are the values for those parameters. To make formatting
easier, a template file is generated for you, just replace the values and add/remove
rows as necessary. See the section for the [!badge corners="pill" text="impute"](/Workflows/impute.md)
module for details on these parameters. The template file will look like:

```text params.stitch
name model	usebx	bxlimit	k	s	ngen
k10_ng50 diploid	TRUE	50000	3	2	10
k1_ng30 diploid	TRUE	50000	3	1	5
high_ngen   diploid TRUE    50000   15  1   100
```
---

## resume
When calling a workflow (e.g. [!badge corners="pill" text="qc"](qc.md)), Harpy performs various file checks
and validations, sets up the Snakemake command, output folder(s), etc. In the event you want to continue a
failed or manually terminated workflow without overwriting the workflow files (e.g. `config.harpy.yaml`),
you can use [!badge corners="pill" text="harpy resume"]. Using `resume` also skips all input/argument validations.

```bash usage
harpy resume [--conda] DIRECTORY
```

#### arguments
{.compact}
| argument    | description                                                                            |
| :---------- | :------------------------------------------------------------------------------------- |
| `DIRECTORY` | [!badge variant="info" text="required"] Output directory of an existing harpy workflow |
| `--conda`   | Generate a `/workflow/envs` folder with the necessary conda enviroments                |

The `DIRECTORY` is the output directory of a previous harpy-invoked workflow, which **must** have the `workflow/config.yaml` file.
For example, if you previously ran `harpy align bwa -o align-bwa ...`, then you would use `harpy resume align-bwa`,
which would have the necessary `workflow/config.yaml` (and other necessary things) required to successfully continue the workflow.
Using [!badge corners="pill" text="resume"] does **not** overwrite any preprocessing files in the target directory (whereas rerunning the workflow would),
which means you can also manually modify the `config.yaml` file (advanced, not recommended unless you are confident with what you're doing).

[!badge corners="pill" text="resume"] also requires an existing and populated `workdir/envs/` directory in the target directory, like the kind all
main `harpy` workflows would create. If one is not present, you can use `--conda` to create one.


---

## view
This convenience command lets you view important details within a Harpy workflow directory
without having to fish around for the right files.

```bash usage
harpy view MODE DIRECTORY
```

```bash example
harpy view config Align/bwa
```

### arguments
{.compact}
| argument      |  description                                                                           |
| :------------ | :------------------------------------------------------------------------------------- |
| `DIRECTORY`   | [!badge variant="info" text="required"] Output directory of an existing harpy workflow |

### modes
{.compact}
| MODE           | description                                                       |
|:---------------|:------------------------------------------------------------------|
| `config`       | View a workflow's config file                                     |
| `environments`* | View the conda environments and their software in `.environments` |
| `log`          | View a workflow's last log file                                   |
| `snakefile`    | View a workflow's snakefile                                       |
| `snakeparams`  | View a workflow's snakemake parameter file                        |

#### environments
`view environments` is an exception in that it does not require a `DIRECTORY` argument.
To use `environments`, you can run it without arguments or give it a `SOFTWARE` argument
to only print environments where `SOFTWARE` is found (it also works with partial matches).

``` viewing all environments
harpy view environments

.environments/f5cb053d77e72fca1c7b6463448fd855_
  - falco=1.2.4
  - fastp
  - multiqc=1.28
  - pysam=0.22
  - quickdeconvolution

.environments/9f9995de2cc77c81654693b8fb002922_
  - quarto
  - r-dt
  - r-dplyr
  - r-highcharter
  - r-magrittr
  - r-plotly
  - r-scales
  - r-tidyr
  - r-viridislite
  - r-xml2
  - r-biocircos
```

``` searching for a program
harpy view environments fast

.environments/f5cb053d77e72fca1c7b6463448fd855_
  - falco=1.2.4
  → fastp
  - multiqc=1.28
  - pysam=0.22
  - quickdeconvolution
```