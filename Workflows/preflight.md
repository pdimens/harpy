---
label: Preflight
description: Run file format checks on haplotagged FASTQ/BAM files
icon: rocket
order: 5
---

# :icon-rocket: Pre-flight checks for input files

===  :icon-checklist: You will need
- at least 2 cores/threads available
- [!badge corners="pill" text="preflight bam"]: SAM/BAM alignment files [!badge variant="secondary" text="BAM recommended"]
- [!badge corners="pill" text="preflight fastq"]: paired-end reads from an Illumina sequencer in FASTQ format [!badge variant="secondary" text="gzip recommended"]
    - **sample names**: [!badge variant="success" text="a-z"] [!badge variant="success" text="0-9"] [!badge variant="success" text="."] [!badge variant="success" text="_"] [!badge variant="success" text="-"] [!badge variant="secondary" text="case insensitive"]
    - **forward**: [!badge variant="success" text="_F"] [!badge variant="success" text=".F"] [!badge variant="success" text=".1"] [!badge variant="success" text="_1"] [!badge variant="success" text="_R1_001"] [!badge variant="success" text=".R1_001"] [!badge variant="success" text="_R1"] [!badge variant="success" text=".R1"] 
    - **reverse**: [!badge variant="success" text="_R"] [!badge variant="success" text=".R"] [!badge variant="success" text=".2"] [!badge variant="success" text="_2"] [!badge variant="success" text="_R2_001"] [!badge variant="success" text=".R2_001"] [!badge variant="success" text="_R2"] [!badge variant="success" text=".R2"] 
    - **fastq extension**: [!badge variant="success" text=".fq"] [!badge variant="success" text=".fastq"] [!badge variant="secondary" text="case insensitive"]
===

Harpy does a lot of stuff with a lot of software and each of these programs expect the incoming data to follow particular formats (plural, unfortunately).
These formatting opinions/specifics are at the mercy of the original developers and while there are times when Harpy can (and does)
modify input/output files for format compatability, it's not always feasible or practical to handle all possible cases. So, our
solution is perform what we lovingly call "pre-flight checks" to assess if your input FASTQ or BAM files are formatted correctly
for the pipeline. There are separate [!badge corners="pill" text="preflight fastq"] and [!badge corners="pill" text="preflight bam"] submodules and the result of each is a report detailing file format quality checks. 

## When to run
- [!badge corners="pill" text="preflight fastq"]: the preflight checks for FASTQ files are best performed _after_ demultiplexing (or trimming/QC) and _before_ sequence alignment
- [!badge corners="pill" text="preflight bam"]: the preflight checks for BAM files should be run _after_ sequence alignment and _before_ consuming those files for other purposes
(e.g. variant calling, phasing, imputation)


```bash fastq usage and example
harpy preflight fastq OPTIONS... INPUTS...

# example 
harpy preflight fastq --threads 20 raw_data
```

```bash bam usage and example
harpy preflight bam OPTIONS... INPUTS...

# example
harpy preflight bam --threads 20 Align/bwa
```

## :icon-terminal: Running Options
In addition to the [!badge variant="info" corners="pill" text="common runtime options"](/commonoptions.md), the [!badge corners="pill" text="preflight fastq"] and [!badge corners="pill" text="preflight bam"] modules are configured using only command-line input arguments:

{.compact}
| argument | short name | description                                                                                                                           |
| :------- | :--------: | :------------------------------------------------------------------------------------------------------------------------------------ |
| `INPUTS` |            | [!badge variant="info" text="required"] Files or directories containing [input fastq or bam files](/commonoptions.md#input-arguments) |

## Workflow

+++ fastq checks
Below is a table of the format specifics [!badge corners="pill" text="preflight fastq"] checks for FASTQ files. Since 10X data doesn't use
the haplotagging data format, you will find little value in running [!badge corners="pill" text="preflight fastq"] on 10X FASTQ files. Take note
of the language such as when "any" and "all" are written.

{.compact}
 Criteria            | Pass Condition                                                               | Fail Condition                                                |
 :------------------ | :--------------------------------------------------------------------------- | :------------------------------------------------------------ |
 AxxCxxBxxDxx format | **all** reads with BX:Z: tag have properly formatted `AxxCxxBxxDxx` barcodes | **any** BX:Z: barcodes have incorrect format                  |
 follows SAM spec    | **all** reads have proper `TAG:TYPE:VALUE` comments                          | **any** reads have incorrectly formatted comments             |
 BX:Z: last comment  | **all** reads have `BX:Z`: as final comment                                  | **at least 1 read** doesn't have `BX:Z:` tag as final comment |
 BX:Z: tag           | any `BX:Z:` tags present                                                     | **all** reads lack `BX:Z:` tag                                |

+++ bam checks
Below is a table of the format specifics [!badge corners="pill" text="preflight bam"] checks for SAM/BAM files. Take note
of the language such as when "any" and "all" are written.

{.compact}
| Criteria            | Pass Condition                                                                    | Fail Condition                                                |
| :------------------ | :-------------------------------------------------------------------------------- | :------------------------------------------------------------ |
| name matches        | the file name matches the `@RG ID:` tag in the header                             | file name does not match `@RG ID:` in the header              |
| MI: tag             | **any** alignments with `BX:Z:` tags also have `MI:i:` (or `MI:Z:`) tags          | **all** reads have `BX:Z:` tag present but `MI:i:` tag absent |
| BX:Z: tag           | any `BX:Z:` tags present                                                          | **all** alignments lack `BX:Z:` tag                           |
| AxxCxxBxxDxx format | **all** alignments with BX:Z: tag have properly formatted `AxxCxxBxxDxx` barcodes | **any** `BX:Z:` barcodes have incorrect format                |
| BX:Z: last tag      | **all** reads have `BX:Z`: as final tag in alignment records                      | **at least 1 read** doesn't have `BX:Z:` tag as final tag     |

+++ output
The default output directory is `Preflight/fastq` or `Preflight/bam` depending on which mode you are using.

+++ Reports
The result of `preflight` is a single HTML report in `inputdir/Preflight/filecheck.xxx.html` where `xxx` is either `fastq` or `bam`
depending on which filetype you specified. The reports for both `fastq` and `bam` are very similar and give you both the
criteria of what type of format checking occurred, the context, relevance, and severity of those checks, along with pass/fails for each
file (or sample).

||| FASTQ file report
![Preflight/filecheck.fastq.html](/static/report_preflightfastq.png)
||| BAM file report
![Preflight/filecheck.bam.html](/static/report_preflightbam.png)
|||
+++
