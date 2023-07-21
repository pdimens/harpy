---
label: Common Options
icon: list-unordered
order: 96
---

# :icon-list-unordered: Common Harpy Options
Every Harpy module has a series of configuration parameters. These are arguments you need to input
to configure the module to run on your data, such as the directory with the reads/alignments,
the genome assembly, etc. All modules (except `extra`) also share a series of common runtime
parameters that don't impact the results of the module, but instead control the speed/verbosity/etc.
of calling the module. These runtime parameters are listed in the modules' help strings and can be 
configured using these arguments:

| argument       | short name | type        |    default    | required | description                                                                                     |
|:---------------|:----------:|:------------|:-------------:|:--------:|:------------------------------------------------------------------------------------------------|
| `--threads`    |    `-t`    | integer     |       4       |    no    | Number of threads to use                                                                        |
| `--snakemake`  |    `-s`    | string      |               |    no    | Additional [Snakemake](snakemake/#adding-snakamake-parameters) options, in quotes |
| `--quiet`      |    `-q`    | toggle      |               |    no    | Supressing Snakemake printing to console                                                        |
| `--help`       |            |             |               |          | Show the module docstring                                                                       |

As as example, you could call the `harpy align` module and specify 20 threads with no output to console:
```bash combining config and runtime arguments
harpy align --threads 20 --directory samples/trimmedreads --method bwa --quiet

# same as #

harpy align -t 20 -d samples/trimmedreads -m bwa -q
```