---
label: Common Options
icon: list-unordered
order: 4
---

# :icon-list-unordered: Common Harpy Options
## Input Arguments
Each of the main Harpy modules (e.g. [!badge corners="pill" text="qc"](Modules/qc.md) or [!badge corners="pill" text="phase"](Modules/phase.md)) follows the format of
```bash
harpy module options arguments
```
where `module` is something like [!badge corners="pill" text="impute"](Modules/impute.md) or [!badge corners="pill" text="snp mpileup"](Modules/snp.md) and `options` are the runtime parameters,
which can include things like an input `--vcf` file, `--molecule-distance`, etc. After the options
is where you provide the input files/directories without flags and following standard BASH expansion
rules (e.g. wildcards). You can mix and match entire directories, individual files, and wildcard expansions.
In most cases, you can provide an unlimited amount of input arguments, which Harpy will parse and symlink
into the `*/workflow/input` folder, leaving the original files unmodified. In practice, that can look like:
```bash
harpy align bwa -t 5 -g genome.fasta data/pop1 data/pop2/trimmed*gz data/pop3/sample{1,2}* data/pop4/sample{2..5}*gz 
```
!!!info not recursive
Keep in mind that Harpy will not recursively scan input directories for files. If you provide `data/` as an input,
Harpy will search for fastq/bam files in `data/` and not in any subdirectories within `data/`. This is done deliberately
to avoid unexpected behavior.
!!!

!!!warning clashing names
Given the regex pattern matching Harpy employs under the hood and the isolation of just the sample names for Snakemake rules,
files in different directories that have the same name (ignoring extensions) will clash. For example, `lane1/sample1.F.fq`
and `lane2/sample1.F.fq` would both derive the sample name `sample1`, which, in a workflow like [!badge corners="pill" text="align"](Modules/Align/Align.md)
would both result in `output/sample1.bam`, creating a problem. This also holds true for the same sample name but different extension, such
as `sample1.F.fq` and `sample1_R1.fq.gz`, which would again derive `sample1` as the sample name and create a naming clash for workflow outputs.
During parsing, Harpy will inform you of naming clashes and terminate to protect you against this behavior. 
!!!

## Common command-line options
Every Harpy module has a series of configuration parameters. These are arguments you need to input
to configure the module to run on your data, such as the directory with the reads/alignments,
the genome assembly, etc. All main modules (e.g. [!badge corners="pill" text="qc"](Modules/qc.md)) also share a series of common runtime
parameters that don't impact the results of the module, but instead control the speed/verbosity/etc.
of calling the module. These runtime parameters are listed in the modules' help strings and can be 
configured using these arguments:

{.compact}
| argument        | short name | type    | default | description                                                                       |
|:--------------- |:----------:|:------- |:-------:|:--------------------------------------------------------------------------------- |
| `--output-dir`  |   `-o`     | string  | varies  | Name of output directory                                                          |
| `--threads`     |   `-t`     | integer | 4       | Number of threads to use                                                          |
| `--conda`       |            | toggle  |         | Use local conda environments instead of preconfigured Singularity container       |
| `--skipreports` |            | toggle  |         | Skip the processing and generation of HTML reports in a workflow                  |
| `--snakemake`   |            | string  |         | Additional [Snakemake](snakemake/#adding-snakamake-parameters) options, in quotes |
| `--quiet`       |   `-q`     | toggle  |         | Suppress Snakemake printing to console                                            |
| `--help`        |   `-h`     |         |         | Show the module docstring                                                         |

As as example, you could call [!badge corners="pill" text="align strobe"](Modules/Align/strobe.md) and specify 20 threads with no output to console:

```bash
harpy align strobe --threads 20 --quiet samples/trimmedreads

# identical to #

harpy align strobe -t 20 -q samples/trimmedreads
```
---

## The `workflow` folder
When you run one of the main Harpy modules, the output directory will contain a `workflow` folder. This folder is
both necessary for the module to run and is very useful to understand what the module did, be it for your own
understanding or as a point of reference when writing the Methods within a manuscript. The presence of the folder
and the contents therein also allow you to rerun the workflow manually. The `workflow` folder may contain the following:

{.compact}
| item | contents | utility |
|:-----|:---------|:--------|
|`*.smk`               | Snakefile with the full recipe of the workflow | useful for understanding the workflow |
| `config.yml`         | Configuration file generated from command-line arguments and consumed by the Snakefile | useful for bookkeeping | 
| `input/`             | Symlinks to all of the provided input files with standardized extensions |
| `report/*.Rmd`       | RMarkdown files used to generate the fancy reports | useful to understand math behind plots/tables or borrow code from |
| `*.summary` | Plain-text overview of the important parts of the workflow | useful for bookkeeping and writing Methods |

---

## The `Genome` folder
You will notice that many of the workflows will create a `Genome` folder in the working 
directory. This folder is to make it easier for Harpy to store the genome and the associated
indexing/etc. files across workflows without having to redo things unnecessarily. Your input 
genome will be symlinked into that directory (not copied, unless a workflow requires gzipping/decompressing),
but all the other files (`.fai`, `.bwt`, `.bed`, etc.) will be created in that directory.
