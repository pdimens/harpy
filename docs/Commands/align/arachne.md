---
label: arachne
description: Align sequences using linked-read information with arachne
category: [linked-read]
tags: [linked-read]
icon: dot
order: 5
---

# :icon-quote: align with linked-read information
Arachne is the successor to lariat, the linked-read aware aligner
originally developed by 10X Genomics. It incorporates linked-read information
to better place alignments, which tends to show improvement when aligning over
repetitive regions. The workflow first preprocesses FASTQ files into the format
expected by Arachne (sorted by barcode, reads with invalid barcodes removed), then
aligns the processed linked-reads with Arachne, while everything else gets aligned
using minibwa.


===  :icon-checklist: You will need
- at least 4 cores/threads available
- a genome assembly in FASTA format: [!badge variant="success" text=".fasta"] [!badge variant="success" text=".fa"] [!badge variant="success" text=".fasta.gz"] [!badge variant="success" text=".fa.gz"] [!badge variant="secondary" text="case insensitive"]
- paired-end fastq sequence files [!badge variant="secondary" icon=":heart:" text="gzipped recommended"]
    - **sample name**: [!badge variant="success" text="a-z"] [!badge variant="success" text="0-9"] [!badge variant="success" text="."] [!badge variant="success" text="_"] [!badge variant="success" text="-"] [!badge variant="secondary" text="case insensitive"]
    - **forward**: [!badge variant="success" text="_F"] [!badge variant="success" text=".F"] [!badge variant="success" text=".1"] [!badge variant="success" text="_1"] [!badge variant="success" text="_R1_001"] [!badge variant="success" text=".R1_001"] [!badge variant="success" text="_R1"] [!badge variant="success" text=".R1"] 
    - **reverse**: [!badge variant="success" text="_R"] [!badge variant="success" text=".R"] [!badge variant="success" text=".2"] [!badge variant="success" text="_2"] [!badge variant="success" text="_R2_001"] [!badge variant="success" text=".R2_001"] [!badge variant="success" text="_R2"] [!badge variant="success" text=".R2"] 
    - **fastq extension**: [!badge variant="success" text=".fq"] [!badge variant="success" text=".fastq"] [!badge variant="secondary" text="case insensitive"]
===

Once sequences have been trimmed and passed through other QC filters, they will need to
be aligned to a reference genome. This module within Harpy expects filtered reads as input,
such as those derived using [!badge corners="pill" text="harpy qc"](../qc.md). You can map reads onto a genome assembly with Harpy using the [!badge corners="pill" text="align arachne"] module:

```bash usage
harpy align arachne OPTIONS... REFERENCE INPUTS...
```
```bash example
harpy align arachne genome.fasta Sequences/ 
```

## :icon-terminal: Running Options
In addition to the [!badge variant="info" corners="pill" text="common runtime options"](/Getting_Started/common_options.md), the [!badge corners="pill" text="align bwa"]/[!badge corners="pill" text="align strobe"] modules are configured using these command-line arguments:

{.compact .clean}
| argument    {.whitespace-nowrap} | default {.whitespace-nowrap} | description                                                                                                                                     |
| :------------------------------- | :--------------------------: | :---------------------------------------------------------------------------------------------------------------------------------------------- |
| `REFERENCE`                      |                              | [!badge variant="info" text="required"] Reference assembly for read mapping                                                                     |
| `INPUTS`                         |                              | [!badge variant="info" text="required"] Files or directories containing [input FASTQ files](/Getting_Started/common_options.md#input-arguments) |
| `-depth-window` `-w`             |           `50000`            | Interval size (in bp) for depth stats                                                                                                           |
| `--extra-params` `-x`            |                              | Additional BWA arguments, in quotes                                                                                                             |
| `--keep-unmapped` `-u`           |            false             | Output unmapped sequences too                                                                                                                   |
| `--molecule-distance` `-d`       |             `0`              | Base-pair distance threshold to separate molecules given as base pairs, disabled with `0`                                                       |
| `--min-quality` `-q`             |             `30`             | Minimum `MQ` (SAM mapping quality) to pass filtering                                                                                            |
