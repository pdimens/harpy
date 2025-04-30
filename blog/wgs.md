---
author: 
  - name: Pavel Dimens
    avatar: https://cals.cornell.edu/sites/default/files/styles/faculty/public/2024-09/afs-headshot-high-res-2cropped_0.jpg
date: 2025-02-06
category: guides
description: How to use Harpy for plain-regular WGS data
icon: feed-rocket
image: https://www.food-safety.com/ext/resources/2021/06/19/WGS_Img01_900.jpg
---

# :icon-feed-rocket: Harpy for (non linked-read) WGS data
As of version `2.0`, Harpy can be used for the early stages of regular whole genome
sequencing (WGS) bioinformatics. Specifically, you can quality checks and trim samples,
align sequences, call SNPs and small indels, phase, and impute genotypes. All of that is done with the flick of
the `--ignore-bx` switch. RADseq data may also work, however the SNP calling workflows
probably won't be very computationally efficient for a highly fragmented RAD assembly. 
There is also another consideration for RADseq regarding marking duplicates (described below).

## Quality Assessment
Using [!badge corners="pill" text="harpy qc"](/Workflows/qc.md), are able to detect and remove adapters, poly G tails, trim low 
quality, bases, detect duplicates with UMIs, etc. You **cannot** use the deconvolution
function of this workflow (`--deconvolve`).
```bash
harpy qc --ignore-bx --trim-adapters auto --min-length 50 data/WGS/sample_*.gz 
```

## Sequence Alignment
Likewise, you can use either [!badge corners="pill" text="harpy align bwa"](/Workflows/Align/bwa.md) or [!badge corners="pill" text="harpy align strobe"](/Workflows/Align/strobe.md) to align
your sequences onto a reference genome. The `--depth-window` and `--molecule-distance`
options are irrelevant and ignored when using `--ignore-bx`. Since EMA is a linked-read
specific aligner, it is not available for WGS/RADseq data, nor would you get any value
from trying to use it for such.

```bash
harpy align bwa --ignore-bx --min-quality 25 genome.fasta data/WGS/trimmed 
```

!!!warning RADseq data
RADseq data will probably work fine too, however you may need to post-process the
BAM files to unset the duplicate flag, as marking duplicates in RADseq (without UMIs) [may cause issues](https://www.researchgate.net/post/How_to_exclude_PCR_duplicates_in_ddRAD) with SNP calling:
```bash
samtools view -b -h --remove-flags 1024 -o output.bam input.bam
```
!!!

## Calling SNPs
The SNP-calling workflows in Harpy don't use linked-read information at all, so you
would use [!badge corners="pill" text="harpy snp mpileup"](/Workflows/snp.md) or [!badge corners="pill" text="harpy snp freebayes"](/Workflows/snp.md) without any modifications.

```bash
harpy snp mpileup --regions 100000 --populations data.groups genome.fasta Align/strobe
```