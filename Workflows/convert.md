---
label: Convert
description: Convert between linked-read data formats
category: [linked-read]
tags: [linked-read]
icon: arrow-switch
nav:
  badge: deprecated|danger
order: 10
---

# :icon-arrow-switch: Convert between data formats [!badge variant="danger" text="being deprecated"]


!!!danger DEPRECATION NOTICE
The convert module is being replaced with [djinn](https://pdimens.github.io/djinn),
which is provided with a conda-based Harpy installation as of version `3.1` (build 2).
It will be officially removed from Harpy starting with version 3.3.
!!!

Regrettably, the bright minds who developed various linked-read technologies cannot seem to agree on a unified data format.
That's annoying at best and hinders the field of linked-read analysis at worst, as there are pieces of very clever software
that are specific to a narrow set of linked-read data formats. Until such a day where there is a concensus, Harpy provides
the means to convert between the various popular linked-read data formats. 

==- Creating order in chaos :warning:
You will notice Harpy offers the "standardize" option, which suggests there is a "standard" format.
That is our attempt to _encourage_ linked-read practioners to consider a unified and practical data format
in which a barcode **of any format** is encoded is in the `BX:Z` tag of a FASTQ/BAM file (a standard SAM-compliant tag) and the validation for the barcode
(whether it is valid or not according to the technology) is encoded in the `VX:i` tag as either `0` (invalid) or `1` (valid).

As an example, if you have stLFR barcoded data, whose barcodes take the form `1_2_3`, the barcode `54_0_1123` would be considered
invalid because stLFR barcodes with a `0` as one of the segments are invalid (missing/ambiguous segment). The `standard` data format,
regardless of FASTQ or BAM, would have the barcoded as `BX:Z:54_0_1123` and the validation as `VX:i:0`.
===

## Convert FASTQ formats
In the event you need your linked-read data converted into a different linked-read format to use outside of Harpy:
don't worry, we got you covered. We might disagree on the fragmented format landscape, but that doesn't mean you
shouldn't be able to use your data how and where you want to. This command converts a paired-end read set of FASTQ
files between the common linked-read FASTQ types.

```bash usage
harpy convert fastq -o <output> -b <barcodes> TARGET FQ1 FQ2
```

```bash example | tellseq → stlfr
harpy convert fastq -o data/orcs_stlfr stlfr data/orcs.R1.fq.gz data/orcs.R2.fq.gz
```

Auto-detects the input format as one of haplotagging, TELLseq, stLFR, or 10X (if `--barcodes` are provided),
and converts it to the format provided as the `TARGET` positional argument. If the input data is
10X format, the `--barcodes` file must contain one nucleotide barcode per line to
determine which barcodes are valid/invalid. In all cases, a file will be created with
the barcode conversion map.

### :icon-move-to-end: Conversion targets

{.compact}
| `TARGET`       | barcode format                                     | example                     |
|:---------------|:---------------------------------------------------|:----------------------------|
| `10x`          | the first N base pairs of R1, given `--barcodes`   |                             |
| `haplotagging` | a `BX:Z:ACBD` SAM tag in the sequence header       | `@SEQID BX:Z:A01C93B56D11`  |
| `stlfr`        | `#1_2_3` format appended to the sequence ID        | `@SEQID#1_2_3`              |
| `tellseq`      | `:ATCG` format appended to the sequence ID         | `@SEQID:GGCAAATATCGAGAAGTC` |


### :icon-terminal: Running Options
{.compact}
| argument          | default | description                                                                   |
|:------------------|:-------:|:------------------------------------------------------------------------------|
| `TARGET`          |         | [!badge variant="info" text="required"] target format for output FASTQ files  |
| `FQ1`             |         | [!badge variant="info" text="required"] forward reads of FASTQ pair           |
| `FQ2`             |         | [!badge variant="info" text="required"] reverse reads of FASTQ pair           |
| `--barcodes` `-b` |         | file of nucleotide barcodes (only necessary when `FROM` is `10x`)             |
| `--output` `-o`   |         | [!badge variant="info" text="required"] file prefix for output fastq files    |
| `--quiet`         |   `0`   | `0` and `1` (all) or `2` (no) output                                          |


## Standardize
In the effort of making it painless to have your data in the preferred standard format, Harpy provides `convert standardize-*`
to quickly standardize FASTQ and BAM files. By default, standardization just moves the barcode (wherever it may be)
into a `BX:Z` SAM tag as-is and does a technology-appropriate validation of the barcode value, which it writes to the
`VX:i` tag. However, you can use `--style` to also convert the barcode style between formats. Keep in mind that each
barcode style has a different upper limit as to how many unique barcodes it can support, which may prevent successful conversions.
The styles are given as:

{.compact}
| Style          | Maximum Unique | What they look like | Example            |
|:---------------|:--------------:|:--------------------|:-------------------|
| `haplotagging` |     $96^4$     | AxxCxxBxxDxx        | A41C22B70D93       |
| `stlfr`        |    $1537^3$    | 1_2_3               | 901_3_1121         |
| `tellseq`      |    $4^{18}$    | 18-base nucleotide  | AGCCATGTACGTATGGTA |
| `10X`          |    $4^{16}$    | 16-base nucleotide  | GGCTGAACACGTGCAG   |

### BAM
If barcodes are present in the sequence name (stlfr, tellseq), this method moves the barcode to the `BX:Z`
tag of the alignment, maintaining the same barcode style by default (auto-detected). If moved to or already in a `BX:Z` tag,
will then write a complementary `VX:i` tag to describe barcode validation `0` (invalid) or `1` (valid).
Use `--style` to also convert the barcode to a different style (`haplotagging`, `stlfr`, `tellseq`, `10X`),
which also writes a `conversion.bc` file to the working directory mapping the barcode conversions. Writes to `stdout`.

```bash usage
harpy convert standardize-bam [--quiet --style] SAM > output.bam
```

```bash example | standardize a bam and change the barcodes to stLFR style
harpy convert standardize-bam --style stflr yucca.bam > yucca.std.stlfr.bam
```

#### :icon-terminal: Running Options
{.compact}
| argument       | default | description                                                                        |
|:---------------|:-------:|:-----------------------------------------------------------------------------------|
| `SAM`          |         | [!badge variant="info" text="required"] input SAM/BAM alignment file               |
| `--quiet`      |   `0`   | `0` and `1` (all) or `2` (no) output                                               |
| `-s`/`--style` |         | change barcode style in the output BAM: [`10x`,`haplotagging`, `stlfr`, `tellseq`] |

### FASTQ
This conversion moves the barcode to the `BX:Z` tag in fastq records, maintaining the same barcode type by default (auto-detected).
See [this section](/Getting_Started/linked_read_data.md#linked-read-data-types) for the location and format expectations for different linked-read technologies.
Also writes a `VX:i` tag to describe barcode validation `0` (invalid) or `1` (valid).
Use `--style` to also convert the barcode to a different style (`haplotagging`, `stlfr`, `tellseq`, `10X`),
which will also write a `conversion.bc` file to the working directory mapping the barcode conversions.

!!!warning Incompatible with 10X data
Standardization will **not** work with the 10X FASTQ format, where the barcodes are the first 16 bases of read 1.
Instead, use [!badge corners="pill" text="convert fastq"](#convert-fastq-formats) to first convert
it into another format like TELLseq. You can standardize it into a barcode style of your liking afterwards.
!!!

```bash usage
harpy convert standardize-fastq [--quiet --style] PREFIX R1.fq R2.fq
```

```bash example | standardize a fastq pair and change the barcodes to stLFR style
harpy convert standardize-fastq --style stflr myotis.stlfr myotis.R1.fq.gz myotis.R2.fq.gz
```

#### :icon-terminal: Running Options
{.compact}
| argument       | default | description                                                                          |
|:---------------|:-------:|:-------------------------------------------------------------------------------------|
| `R1.fq`        |         | [!badge variant="info" text="required"] input FASTQ forward-read file                |
| `R1.fq`        |         | [!badge variant="info" text="required"] input FASTQ reverse-read file                |
| `PREFIX`       |         | [!badge variant="info" text="required"] prefix for output filenames                  |
| `--quiet`      |   `0`   | `0` and `1` (all) or `2` (no) output                                                 |
| `-s`/`--style` |         | change barcode style in the output FASTQ: [`10x`,`haplotagging`, `stlfr`, `tellseq`] |


----

:::info Useless trivia
The original version of this command was written while Pavel was waiting at a mechanic shop for his car to be repaired. During development,
it was called `lr-switcheroo`.
:::