---
label: Common Options
icon: list-unordered
order: 5
---

# :icon-list-unordered: Common Harpy Options
## Input Arguments
Each of the main Harpy modules (e.g. [!badge corners="pill" text="qc"](/Workflows/qc.md) or [!badge corners="pill" text="phase"](/Workflows/phase.md)) follows the format of
```bash
harpy command options arguments
```
where `command` is something like [!badge corners="pill" text="impute"](/Workflows/impute.md) or [!badge corners="pill" text="snp mpileup"](/Workflows/snp.md) and `options` are the runtime parameters,
which can include things like `--molecule-distance`, etc. After the options
is where you provide the input files/directories without flags and following standard BASH expansion
rules (e.g. wildcards). You can mix and match entire directories, individual files, and wildcard expansions.
In most cases, you can provide an unlimited amount of input arguments. In practice, that can look like:
```bash
harpy align bwa -t 5 genome.fasta data/pop1 data/pop2/trimmed*gz data/pop3/sample{1,2}* data/pop4/sample{2..5}*gz 
```
!!!info not recursive
By design, Harpy will not recursively scan input directories for files. If you provide `data/` as an input,
Harpy will search for fastq/bam files in `data/` and not in any subdirectories within `data/`. This is done
to avoid unexpected behavior.
!!!

!!!warning clashing names
Given the regex pattern matching that happens under the hood and the isolation of just the sample names for Snakemake rules,
files in different directories that have the same name (ignoring extensions) will clash. For example, `lane1/sample1.F.fq`
and `lane2/sample1.F.fq` would both derive the sample name `sample1`, which, in a workflow like [!badge corners="pill" text="align"](/Workflows/Align/Align.md)
would both result in `output/sample1.bam`, creating a problem. This also holds true for the same sample name but different extension, such
as `sample1.F.fq` and `sample1_R1.fq.gz`, which would again derive `sample1` as the sample name and create a naming clash for workflow outputs.
You will be informed ahead of time of and naming clashes to protect you against this happening. 
!!!

## Software Dependencies
Harpy workflows typically require various different pieces of software to run. To
keep the Harpy installation small, we include only the bare minimum to invoke Harpy.
Everything else (e.g. `freebayes`, `hapcut2`, etc.) is installed as needed at runtime by Snakemake.
By default, Harpy has Snakemake to install a workflow's software dependencies as local conda environments
in the `.environments` folder, however you can use `--container` to instead have Snakemake use a pre-configured
Harpy container hosted on Dockerhub to manage workflow dependencies.

## Common command-line options
Every Harpy module has a series of configuration parameters. These are arguments you need to input
to configure the module to run on your data, such as the directory with the reads/alignments,
the genome assembly, etc. All main modules (e.g. [!badge corners="pill" text="qc"](/Workflows/qc.md)) also share a series of common runtime
parameters that don't impact the results of the module, but instead control the speed/verbosity/etc.
of calling the module. These runtime parameters are listed in the modules' help strings and can be 
configured using these arguments:

{.compact}
| argument            | type              | default | description                                                                                                                                            |
|:--------------------|:------------------|:-------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--container`       | toggle            |         | Use preconfigured Apptainer container instead of local conda environments                                                                            |
| `--contigs`         | file path or list |         | Contigs to plot in the report(s)                                                                                                                       |
| `--help` `-h`       |                   |         | Show the module docstring                                                                                                                              |
| `--hpc`             |                   |         | Have snakemake submit all jobs to an HPC ([details](Resources/hpc.md))                                                                                 |
| `--output-dir` `-o` | string            | varies  | Name of output directory                                                                                                                               |
| `--quiet`           | integer [0,1,2]   |    0    | `0` prints all progress information, `1` prints unified progress bar, `2` suppressess all console output except errors                                 |
| `--setup-only`      | toggle            |         | [!badge variant="secondary" corners="pill" text="hidden"](/Workflows/qc.md) Perform validations and setup workflow environment, but don't run anything |
| `--skip-reports`    | toggle            |         | Skip the processing and generation of HTML reports in a workflow                                                                                       |
| `--snakemake`       | string            |         | Additional [Snakemake](Resources/snakemake#adding-snakemake-parameters) options, in quotes                                                             |
| `--threads` `-t`    | integer           |    4    | Number of threads to use                                                                                                                               |
| `--unlinked` `-U`   | toggle            |         | Treat the input as non linked-read data                                                                                                                |

### --unlinked
As of version 3.0, Harpy tries to auto-detect the type of data that input FASTQ or BAM files may be: `haplotagging`, `stlfr`, `tellseq`, or `none`.
The formats for these data and how their barcodes look are described [in this section](linked_read_data.md#linked-read-data-types). That means
you no longer have to specify linked-read technology, and if you aren't using linked reads, that's fine too ([in most cases](Guides/wgs_data.md)). 
However, you can force many commands to treat your data as not linked-read using `-U`/`--unlinked`. This toggle skips data type detection
and circumvents any linked-read specific parts of workflows. Data that is already not linked-read would be detected as such, so this
tends to be more helpful if you want to either 1. shave off a second or two from preprocessing non-linked-read data or 2. push linked-read
data through a workflow pretending it's not linked-read data.

### --contigs
Some of the workflows (like [!badge corners="pill" text="align"](/Workflows/Align/Align.md)) plot per-contig information in their reports.
By default, Harpy will plot up to 30 of the largest contigs. If you are only interested in a specific set of contigs, then you can use `--contigs`
to have Harpy only create plots for those contigs. **This will only impact plotting for reports**. This can be done by including a file of one-per-line contig names or a comma-separated
list of contigs (without spaces):

==- contigs.txt
```
contig1
contig2
sexchrom1
```
===
```bash
harpy align bwa --contigs contig1,contig2,sexchrom1 genome.fasta dir/data
# or #
harpy align bwa --contigs contigs.txt genome.fasta dir/data
```
!!!warning too many contigs
Things start to look sloppy and cluttered with >30 contigs, so it's advisable not to
exceed that number.
!!!

### example
You could call [!badge corners="pill" text="align strobe"](/Workflows/Align/strobe.md) and specify 20 threads with no output to console:

```bash
harpy align strobe --threads 20 --quiet 2 genome.fasta samples/trimmedreads

# identical to #

harpy align strobe -t 20 -q genome.fasta samples/trimmedreads
```
---

## The `workflow` folder
When you run one of the main Harpy modules, the output directory will contain a `workflow` folder. This folder is
both necessary for the module to run and is very useful to understand what the module did, be it for your own
understanding or as a point of reference when writing the Methods within a manuscript. The presence of the folder
and the contents therein also allow you to rerun the workflow manually. The `workflow` folder may contain the following:

{.compact}
| item               | contents                                                                                                       | utility                                                |
|:-------------------|:---------------------------------------------------------------------------------------------------------------|:-------------------------------------------------------|
| `*.smk`            | Snakefile with the full recipe of the workflow                                                                 | understanding the entire workflow                      |
| `config.yml`       | Configuration file for Snakemake workflow dispatching                                                          | general bookkeeping, advanced runs                     |
| `config.harpy.yml` | Configuration file generated from command-line arguments and consumed by the Snakefile                         | general bookkeeping, advanced runs                     |
| `envs/`            | Configurations of the software environments required by the workflow                                           | bookkeeping                                            |
| `hpc/`             | Folder with the HPC-specific configuration file that let's Snakemake submit jobs to a scheduler on your behalf | necessary for running on an HPC                        |
| `reference/`       | Folder with a link or copy to the FASTA file used as the reference for various workflows                       | necessary for concurrent workflows to avoid data races |
| `report/*.qmd`     | Quarto files used to generate the fancy reports                                                                | seeing math behind plots/tables or borrow code from    |
| `*.summary`        | Plain-text overview of the important parts of the workflow                                                     | bookkeeping and writing Methods in manuscripts         |
