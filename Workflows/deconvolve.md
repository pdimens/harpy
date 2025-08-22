---
label: Deconvolve
description: Resolve barcodes shared by different molecules
category: [linked-read]
tags: [linked-read]
icon: tag
order: 10
---

# :icon-tag: Resolve barcodes shared by different molecules

===  :icon-checklist: You will need
- paired-end reads from an Illumina sequencer in FASTQ format [!badge variant="secondary" icon=":heart:" text="gzipped recommended"]
    - **sample name**: [!badge variant="success" text="a-z"] [!badge variant="success" text="0-9"] [!badge variant="success" text="."] [!badge variant="success" text="_"] [!badge variant="success" text="-"] [!badge variant="secondary" text="case insensitive"]
    - **forward**: [!badge variant="success" text="_F"] [!badge variant="success" text=".F"] [!badge variant="success" text=".1"] [!badge variant="success" text="_1"] [!badge variant="success" text="_R1_001"] [!badge variant="success" text=".R1_001"] [!badge variant="success" text="_R1"] [!badge variant="success" text=".R1"] 
    - **reverse**: [!badge variant="success" text="_R"] [!badge variant="success" text=".R"] [!badge variant="success" text=".2"] [!badge variant="success" text="_2"] [!badge variant="success" text="_R2_001"] [!badge variant="success" text=".R2_001"] [!badge variant="success" text="_R2"] [!badge variant="success" text=".R2"] 
    - **fastq extension**: [!badge variant="success" text=".fq"] [!badge variant="success" text=".fastq"] [!badge variant="secondary" text="case insensitive"]
===



Running [!badge corners="pill" text="deconvolve"] is **optional**. In the alignment
workflows ([!badge corners="pill" text="align bwa"](Align/bwa.md) 
[!badge corners="pill" text="align strobe"](Align/strobe.md)), Harpy already uses a distance-based approach to
deconvolve barcodes and assign `MI` tags (Molecular Identifier). This workflow uses a reference-free method,
[QuickDeconvolution](https://github.com/RolandFaure/QuickDeconvolution), which uses k-mers to look at "read clouds" (all reads with the same linked-read barcode)
and decide which ones likely originate from different molecules. Regardless of whether you run 
this workflow or not, [!badge corners="pill" text="harpy align"](Align/Align.md) will still perform its own deconvolution.

!!! Also in harpy qc
This method of deconvolution is also available as an option in the [!badge corners="pill" text="qc"](qc.md) workflow
!!!

```bash usage
harpy deconvolve OPTIONS... INPUTS...
```

```bash example | deconvolve with default parameters
harpy deconvolve path/to/data/*.fq
```

## :icon-terminal: Running Options
{.compact}
| argument             | default | description                                                                                                                     |
|:---------------------|:-------:|:--------------------------------------------------------------------------------------------------------------------------------|
| `INPUTS`             |         | [!badge variant="info" text="required"] Files or directories containing [input FASTQ files](/Getting_Started/common_options.md#input-arguments) |
| `--density` `-d`     |   `3`   | On average, $\frac{1}{2^d}$ kmers are indexed                                                                                   |
| `--dropout` `-a`     |   `0`   | Minimum cloud size to deconvolve                                                                                                |
| `--kmer-length` `-k` |  `21`   | Size of k-mers to search for similarities                                                                                       |
| `--window-size` `-w` |  `40`   | Size of window guaranteed to contain at least one kmer                                                                          |

## Resulting Barcodes
After deconvolution, some barcodes may have a hyphenated suffix like `-1` or `-2` (e.g. `A01C33B41D93-1`).
This is how deconvolution methods create unique variants of barcodes to denote that identical barcodes
do not come from the same original molecules. QuickDeconvolution adds the `-0` suffix to barcodes it was unable
to deconvolve.

## Harpy Deconvolution Nuances
Some of the downstream linked-read tools Harpy uses expect linked read barcodes to either look like the 16-base 10X
variety or a standard haplotag (AxxCxxBxxDxx). Their pattern-matching would not recognize barcodes deconvoluted with
hyphens. To remedy this, `MI` assignment in [!badge corners="pill" text="align bwa"](Align/bwa.md)
and [!badge corners="pill" text="align strobe"](Align/strobe.md) will assign the deconvolved (hyphenated) barcode to a `DX:Z`
tag and restore the original barcode as the `BX:Z` tag.
