---
label: Simulate
description: Simulate genomic data
---

# :icon-flame: Simulate Genomic Data

You may be interested in benchmarking variant detection or maybe just trying out
haplotagging data without any financial commitment-- that's where simulations
come in handy. 

## Simulate Genomic Variants
Harpy lets you simulate genomic variants via [!badge corners="pill" text="simulate {snpindel,inversion,...}"](simulate-variants.md) for different variant
types such as single nucleotide polymorphisms (SNP), indels, inversions, copy number variants (CNV), and translocations. All you need is to provide a genome to simulate
those variants onto.

## Simulate Haplotag Linked-Reads
You can also simulate haplotag-style linked reads from an existing genome using [!badge corners="pill" text="simulate linkedreads"](simulate-linkedreads.md). Harpy
incorporates [LRSIM](https://github.com/aquaskyline/LRSIM) to generate linked reads
from a diploid genomic. If you only have a haploid genome, then you can create a diploid genome by simulating variants into it with [!badge corners="pill" text="simulate {snpindel,inversion,...}"](simulate-variants.md).