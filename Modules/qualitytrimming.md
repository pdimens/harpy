---
label: Quality Trimming
icon: codescan-checkmark
order: 6
---

# Quality Trimming Sequence Data
|||  :icon-checklist: You will need
- at least 2 cores/threads available
- paired-end b/gzipped fastq sequence files
|||

You can remove adapters and quality trim sequences using:
```bash
harpy trim OPTIONS... 
```

## Running Options
| argument         | short name | type        | default | required | description                                                            |
|:-----------------|:----------:|:------------|:-------:|:--------:|:-----------------------------------------------------------------------|
| `--dir`          |    `-d`    | folder path |         | **yes**  | Directory with sequence alignments                                     |
| `--max-length`   |    `-l`    | integer     |   150   |    no    | Maximum length to trim sequences down to                               |
| `--extra-params` |    `-x`    | string      |         |    no    | Additional Hapcut2 parameters, in quotes                               |
| `--threads`      |    `-t`    | integer     |    4    |    no    | Number of threads to use                                               |
| `--snakemake`    |    `-s`    | string      |         |    no    | Additional Snakemake options, in quotes ([more info](../snakemake.md)) |
| `--help`         |            |             |         |          | Show the module docstring                                              |

## Fastq file format
There are a handful of "accepted" naming schemes for fastq file extensions, but Harpy only accepts a limited number of them, shown below.
The fastq files **must** be bzipped or gzipped and be **consistent** with regards to the extensions and read-pair naming styles.
That is, all your files must only use `.fastq.gz` or only use `.fq.gz` for all files, and the same for `.1.`/`.2.` or `.F.`/`.R.`.
Notice that the read pair part differs from the [accepted fastq formats](readmapping.md/#fastq-file-format) for aligning reads.
#### acceptable formats
- file extension is either `.fastq.gz` or `.fq.gz` (do not mix)
- forward-reverse is provided as either `.1.`/`.2.` **or** `.F.`/`.R.` (do no mix)
    - _e.g._ `samplename.F.fq.gz` and `samplename.R.fq.gz`
    - _e.g._ `samplename.1.fq.gz` and `samplename.2.fq.gz`
    - or the same but ending with `.fastq.gz`, but don't mix and match

---
## Fastp Workflow
[Fastp](https://github.com/OpenGene/fastp) is an ultra-fast all-in-one adapter remover, deduplicator, 
and quality trimmer. Harpy uses it to remove adapters, low-quality bases, and trim sequences down to a particular
length (default 150bp). Harpy uses the fastp overlap analysis to identify adapters for removal and a sliding window
approach (`--cut-right`) to identify low quality bases. The workflow is quite simple.

```mermaid
graph LR
    A([trim reads]) --> B([create reports])
```