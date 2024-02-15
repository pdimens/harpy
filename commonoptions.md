---
label: Common Options
icon: list-unordered
order: 4
---

# :icon-list-unordered: Common Harpy Options
## Common command-line options
Every Harpy module has a series of configuration parameters. These are arguments you need to input
to configure the module to run on your data, such as the directory with the reads/alignments,
the genome assembly, etc. All main modules (e.g. `qc`) also share a series of common runtime
parameters that don't impact the results of the module, but instead control the speed/verbosity/etc.
of calling the module. These runtime parameters are listed in the modules' help strings and can be 
configured using these arguments:

| argument      | short name | type    | default | required | description                                                                       |
|:------------- |:----------:|:------- |:-------:|:--------:|:--------------------------------------------------------------------------------- |
| `--threads`   | `-t`       | integer | 4       | no       | Number of threads to use                                                          |
| `--print-only` |           | toggle  |         | no       | Perform internal validations, build the workflow environment, and print the resulting Snakemake command, but don't run anything |
| `--skipreports` | `-r`     | toggle  |         | no       | Skip the processing and generation of HTML reports in a workflow                  |
| `--snakemake` | `-s`       | string  |         | no       | Additional [Snakemake](snakemake/#adding-snakamake-parameters) options, in quotes |
| `--quiet`     | `-q`       | toggle  |         | no       | Supressing Snakemake printing to console                                          |
| `--help`      |            |         |         |          | Show the module docstring                                                         |

As as example, you could call the `harpy align` module and specify 20 threads with no output to console:

```bash
harpy align bwa --threads 20 --directory samples/trimmedreads --quiet

# same as #

harpy align bwa -t 20 -d samples/trimmedreads -q
```

## The `workflow` folder
When you run one of the main Harpy modules, the output directory will contain a `workflow` folder. This folder is
both necessary for the module to run and is very useful to understand what the module did, be it for your own
understanding or as a point of reference when writing the Methods within a manuscript. The presence of the folder
and the contents therein also allow you to rerun the workflow manually. The `workflow` folder may contain the following:

| item | contents | utility |
|:-----|:---------|:--------|
|`*.smk`| Snakefile with the full recipe of the workflow | useful for understanding the workflow |
| `config.yml` | Configuration file generated from command-line arguments and consumed by the Snakefile | useful for bookkeeping | 
| `report/*.Rmd` | RMarkdown files used to generate the fancy reports | useful to understand math behind plots/tables or borrow code from |
| `*.workflow.summary` | Plain-text overview of the important parts of the workflow | useful for bookkeeping and writing Methods |

## The `Genome` folder

You will notice that many of the workflows will create a `Genome` folder in the working 
directory. This folder is to make it easier for Harpy to store the genome and the associated
indexing/etc. files. Your input genome will be symlinked into that directory (not copied), but
all the other files (`.fai`, `.bwt`, `.bed`, etc.) will be created in that directory.
