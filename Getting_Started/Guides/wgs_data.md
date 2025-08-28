---
author: 
  - name: Pavel Dimens
    avatar: https://cals.cornell.edu/sites/default/files/styles/faculty/public/2024-09/afs-headshot-high-res-2cropped_0.jpg
date: 2025-02-06
category: guides
description: How to use Harpy for plain-regular WGS data
icon: feed-rocket
label: Non linked-read data
image: https://www.food-safety.com/ext/resources/2021/06/19/WGS_Img01_900.jpg
---

# :icon-feed-rocket: Harpy for (non linked-read) WGS data
As of Harpy v3, the program will auto-detect that your input FASTQ or BAM files are not linked-read data. This can also be forced with `--unlinked` / `-U`.

||| Harpy v2
- When available, use `--lr-type none` to ignore linked-read specific things
- This option was named `--ignore-bx` in versions <2.7

As of version `2.0`, Harpy can be used to process regular whole genome
sequencing (WGS) data. Specifically, you can quality checks and trim samples,
align sequences, call SNPs and small indels, phase, and impute genotypes. All of that is done setting
`--lr-type none` for workflows where `--lr-type` is available (`--ignore-bx` toggle in versions <2.7).
RADseq data may also work, however the SNP calling workflows
probably won't be very computationally efficient for a highly fragmented RAD assembly. 
There is also another consideration for RADseq regarding marking duplicates (described below).

Given some of the reports Harpy produces from its workflows, you can safely ignore the stuff
specific to linked-read information.

## Quality Assessment
Using [!badge corners="pill" text="harpy qc"](/Workflows/qc.md), you are able to detect and remove adapters, poly G tails, trim low 
quality, bases, detect duplicates with UMIs, etc. You **cannot** use `--deconvolve` when ignoring
linked-read information.
```bash qc example
harpy qc --lr-type none --trim-adapters auto --min-length 50 data/WGS/sample_*.gz 
```

## Sequence Alignment
Likewise, you can use either [!badge corners="pill" text="harpy align bwa"](/Workflows/Align/bwa.md) or [!badge corners="pill" text="harpy align strobe"](/Workflows/Align/strobe.md) to align
your sequences onto a reference genome. The `--molecule-distance` will be ignored when
using `--lr-type none`.

```bash align example
harpy align bwa --lr-type none --min-quality 25 genome.fasta data/WGS/trimmed 
```

!!!warning RADseq data
RADseq data will probably work fine too, however you may need to post-process the
BAM files to unset the duplicate flag, as marking duplicates in RADseq (without UMIs) [may cause issues](https://www.researchgate.net/post/How_to_exclude_PCR_duplicates_in_ddRAD) with SNP calling:
```bash clear the duplicate tag
samtools view -b -h --remove-flags 1024 -o output.bam input.bam
```
!!!

## Calling SNPs
The SNP-calling workflows in Harpy don't use linked-read information at all, so you
would use [!badge corners="pill" text="harpy snp mpileup"](/Workflows/snp.md) or [!badge corners="pill" text="harpy snp freebayes"](/Workflows/snp.md) without any modifications.

```bash snp example
harpy snp mpileup --regions 100000 --populations data.groups genome.fasta Align/strobe
```

## Impute Genotypes
You can use the third (`usebx`) column of the [parameter file](/Workflows/impute.md/#parameter-file) to disable the barcode-aware
routines of [!badge corners="pill" text="harpy impute"](/Workflows/impute.md) by setting the value to `FALSE`:
``` stitch.parameters
name    model   usebx   bxlimit   k       s       nGen
model1    diploid   FALSE    50000    10      5       50
model2    diploid   FASE    50000   15      10      100
```
Naturally, ignoring barcodes will also ignore whatever values are set for `bxlimit`. Otherwise, invoke the imputation workflow as you would normally:
```bash impute example
harpy impute -t 10 stitch.parameters data/variants.bcf data/*.bam
```

## Phase Genotypes
Like most of the other workflows, use `--lr-type none` with [!badge corners="pill" text="harpy phase"](/Workflows/phase.md) to perform phasing without incorporating linked-read barcode
information. When using this option, the value for `-d`/`--molecule-distance` will be ignored:
```bash phase example
harpy phase -t 10 --lr-type none variants.bcf data/*.bam 
```

|||
