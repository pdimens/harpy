---
label: Convert
description: Convert between linked-read data formats
category: [linked-read]
tags: [linked-read]
icon: arrow-switch
order: 10
---

# :icon-arrow-switch: Convert between data formats

Regrettably, the bright minds who developed various linked-read technologies cannot seem to agree on a unified data format.
That's annoying at best and hinders the field of linked-read analysis at worst, as there are pieces of very clever software
that are specific to a narrow set of linked-read data formats. Until such a day where there is a concensus, Harpy provides
the means to convert between the various popular linked-read data formats. 

==- Creating order in chaos :warning:
You will notice one of the formats is called
`standard`, and that is our attempt to _encourage_ linked-read practioners to consider a unified and practical data format
in which a barcode **of any format** is encoded is in the `BX:Z` tag of a FASTQ/BAM file (a standard SAM-compliant tag) and the validation for the barcode
(whether it is valid or not according to the technology) is encoded in the `VX:i` tag as either `0` (invalid) or `1` (valid).

As an example, if you have stLFR barcoded data, whose barcodes take the form `1_2_3`, the barcode `54_0_1123` would be considered
invalid because stLFR barcodes with a `0` as one of the segments are invalid (missing/ambiguous segment). The `standard` data format,
regardless of FASTQ or BAM, would have the barcoded as `BX:Z:54_0_1123` and the validation as `VX:i:0`.
===

## Convert FASTQ formats
```bash usage
harpy convert fastq -o <output> -b <barcodes> FROM TO FQ1 FQ2
```

```bash example (tellseq to stlfr)
harpy convert fastq -o data/orcs_tellseq tellseq stlfr data/orcs.R1.fq.gz data/orcs.R2.fq.gz
```

Takes the positional arguments `FROM` to indicate input data format and `TO` is the
target data format. Both of these arguments allow the formats provided in the table below. `10x`
input data requires a `--barcodes` file containing one nucleotide barcode per line to
determine which barcodes are valid/invalid. In all cases, a file will be created with
the barcode conversion map. Requires 2 threads.

{.compact}
| `FROM`/`TO`      | barcode format                                     | example                     |
|:-------------|:---------------------------------------------------|:----------------------------|
| `10x`          | the first N base pairs of R1, given `--barcodes`   |                             |
| `haplotagging` | a `BX:Z:ACBD` SAM tag in the sequence header       | `@SEQID BX:Z:A01C93B56D11`  |
| `standard`     | a `BX:Z` SAM tag in the sequence header, any style | `@SEQID BX:Z:ATAGCAC_AGGA`  |
| `stlfr`        | `#1_2_3` format appended to the sequence ID        | `@SEQID#1_2_3`              |
| `tellseq`      | `:ATCG` format appended to the sequence ID         | `@SEQID:GGCAAATATCGAGAAGTC` |


### :icon-terminal: Running Options
{.compact}
| argument           | default | description                                                                     |
|:-------------------|:-------:|:--------------------------------------------------------------------------------|
| `FROM`             |         | [!badge variant="info" text="required"] current format of the input FASTQ files |
| `TO`               |         | [!badge variant="info" text="required"] target format for output FASTQ files:   |
| `FQ1`              |         | [!badge variant="info" text="required"] forward reads of FASTQ pair             |
| `FQ2`              |         | [!badge variant="info" text="required"] reverse reads of FASTQ pair             |
| `--barcodes` `-b` |         | file of nucleotide barcodes (only necessary when `FROM` is `10x`)               |
| `--output` `-o`   |         | [!badge variant="info" text="required"] file prefix for output fastq files      |
| `--quiet`          |   `0`   | `0` and `1` (all) or `2` (no) output                                            |


## Convert BAM formats
This function converts between linked-read barcode formats in alignments, that is, it
changes the barcode type of the alignment file (SAM/BAM), expecting the barcode to be
in the `BX:Z` tag of the alignment. The barcode type is automatically detected and the
resulting barcode will be in the `BX:Z` tag. Use `--standardize` to optionally standardize
the output file (recommended), meaning a `VX:i` tag is added to describe
barcode validation with `0` (invalid) and `1` (valid). Writes to `stdout`.

{.compact}
| to           | barcode format | example                   |
|:-------------|:---------------|:--------------------------|
| 10x          | nucleotides    | `BX:Z:GGCAAATATCGAGAAGTC` |
| haplotagging | AxxCxxBxxDxx   | `BX:Z:A01C44B81D35`       |
| stlfr        | `1_2_3`        | `BX:Z:911_327_27`         |
| tellseq      | nucleotides    | `BX:Z:GGCAAATATCGAGAAGTC` |

```bash usage
harpy convert bam [--standardize] TO SAM > output.bam
```

```bash example (tellseq to stlfr)
harpy convert bam --standardize haplotagging pomegranate.tellseq.bam > pomegranate.haptag.bam
```

### :icon-terminal: Running Options
{.compact}
| argument        | default | description                                                                                                            |
|:----------------|:-------:|:-----------------------------------------------------------------------------------------------------------------------|
| `TO`            |         | [!badge variant="info" text="required"] barcode format for output BAM file: [`10x`,`haplotagging`, `stlfr`, `tellseq`] |
| `SAM`           |         | [!badge variant="info" text="required"] input SAM/BAM alignment file                                                   |
| `--standardize` | `False` | whether to standardize the output BAM with a `VX:i` validation tag                                                     |
| `--quiet`       |   `0`   | `0` and `1` (all) or `2` (no) output                                                                                   |

## Standardize barcodes
This conversion moves the barcode from the sequence name into the `BX:Z` tag of the alignment,
maintaining the same barcode type (i.e. there is no format conversion). It is intended
for tellseq and stlfr data, which encode the barcode in the read name. Also writes a `VX:i` tag
to describe barcode validation `0` (invalid) or `1` (valid). Writes to `stdout`.

```bash usage
harpy convert standardize [--quiet] SAM > output.bam
```

```bash example
harpy convert standardize yucca.bam > yucca.std.bam
```

### :icon-terminal: Running Options
{.compact}
| argument        | default | description                                                                                                            |
|:----------------|:-------:|:-----------------------------------------------------------------------------------------------------------------------|
| `SAM`           |         | [!badge variant="info" text="required"] input SAM/BAM alignment file                                                   |
| `--quiet`       |   `0`   | `0` and `1` (all) or `2` (no) output                                                                                   |


:::info Useless trivia
This module was written while Pavel was waiting at a mechanic shop for his car to be repaired. During development,
it was called `lr-switcheroo`.
:::