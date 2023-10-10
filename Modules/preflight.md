---
label: Preflight
description: Run file format checks on haplotagged FASTQ/BAM files
icon: rocket
#visibility: hidden
order: 6
---

# :icon-rocket: Pre-flight checks for input files

===  :icon-checklist: You will need
- at least 2 cores/threads available
- `fastq`: paired-end reads from an Illumina sequencer (gzipped recommended)
- `bam`: SAM/BAM alignment files
===

Harpy does a lot of stuff with a lot of software and each of these programs expect the incoming data to follow particular formats (plural, unfortunately).
These formatting opinions/specifics are at the mercy of the original developers and while there are times when Harpy can (and does)
modify input/output files for format compatability, it's not always feasible or practical to handle all possible cases. So, our
solution is perform what we lovingly call "pre-flight checks" to assess if your input FASTQ or BAM files are formatted correctly
for the pipeline. There are separate `fastq` and `bam` submodules and the result of each is a report detailing "file format QC." 

#### when to run
- **BAM**: the preflight checks for BAM files should be run _after_ sequence alignment and _before_ consuming those files for other purposes
(e.g. variant calling, phasing, imputation)

- **FASTQ**: the preflight checks for FASTQ files are best performed _after_ demultiplexing (or trimming/QC) and _before_ sequence alignment

```bash fastq usage and example
harpy preflight fastq OPTIONS...

# example 
harpy preflight fastq --threads 20 -d raw_data
```

```bash bam usage and example
harpy preflight bam OPTIONS... 

# example
harpy preflight bam --threads 20 -d Align/bwa
```

## :icon-terminal: Running Options
In addition to the [common runtime options](../commonoptions.md), the `harpy demultiplex` module is configured using these command-line arguments:

| argument          | short name | type       | default | required | description                                                                          |
|:------------------|:----------:|:-----------|:-------:|:--------:|:-------------------------------------------------------------------------------------|
| `--directory`          |    `-d`    | folder path |         | **yes**  | Directory with sequences or alignments                                                              |

---

## `fastq` checks
Below is a table of the format specifics `harpy preflight` checks for FASTQ files. Take note
of the language such as when "any" and "all" are written.

| Criteria | Pass Condition | Fail Condition |
|:---|:---|:---|
|AxxCxxBxxDxx format| **all** reads with BX:Z: tag have properly formatted `AxxCxxBxxDxx` barcodes | **any** BX:Z: barcodes have incorrect format|
|follows SAM spec | **all** reads have proper `TAG:TYPE:VALUE` comments | **any** reads have incorrectly formatted comments|
|BX:Z: last comment | **all** reads have `BX:Z`: as final comment| **at least 1 read** doesn't have `BX:Z:` tag as final comment|
|BX:Z: tag | any `BX:Z:` tags present | **all** reads lack `BX:Z:` tag|

---

## `bam` checks
Below is a table of the format specifics `harpy preflight` checks for SAM/BAM files. Take note
of the language such as when "any" and "all" are written.

| Criteria | Pass Condition | Fail Condition |
|:---|:---|:---|
|name matches| the file name matches the `@RG ID:` tag in the header| file name does not match `@RG ID:` in the header|
|MI: tag| **any** alignments with `BX:Z:` tags also have `MI:i:` (or `MI:Z:`) tags| **all** reads have `BX:Z:` tag present but `MI:i:` tag absent|
|BX:Z: tag| any `BX:Z:` tags present| **all** alignments lack `BX:Z:` tag|
|AxxCxxBxxDxx format| **all** alignments with BX:Z: tag have properly formatted `AxxCxxBxxDxx` barcodes| **any** `BX:Z:` barcodes have incorrect format|
|BX:Z: last tag| **all** reads have `BX:Z`: as final tag in alignment records | **at least 1 read** doesn't have `BX:Z:` tag as final tag|

## output
Unlike the other modules. `preflight` will not create a new folder in your working directory. Instead, it will create 
a `Preflight` folder in the same directory that was provided for `-d` (`--directory`). This design is intended to keep
the reports near the source data.