---
author: 
  - name: Pavel Dimens
    avatar: https://cals.cornell.edu/sites/default/files/styles/faculty/public/2024-09/afs-headshot-high-res-2cropped_0.jpg
date: 2025-03-18
category: guides
description: What you ought to know about linked reads
icon: git-merge-queue
image: /static/linked_reads.svg
---

# :icon-git-compare: An introduction to linked-read data
Harpy was originally tailor-made for haplotagging linked read data, but now it works for most linked-read technologies
along with non-linked read data. BRL/LRTK/LongRanger are similar pieces of software for Tellseq, stLFR, and 10X linked-read
data. But, what if you don't use linked reads (yet) and want to understand what it actually is? This post walks you through what 
linked-read data is and some of the concepts that are unique to it that makes it different than your typical short-read data.

===  :icon-checklist: You will need
- a cup of tea, coffee, or water
- a can-do attitude
===

## Linked reads
Let's start with the most obvious, what _are_ linked reads, which are sometimes called "synthetic long reads"? They are short read
(e.g. Illumina) data. What makes them different is that they contain an added DNA segment ("barcode") that lets us associate
sequences as having originated from **a single DNA molecule**. That means if we have 4 sequences that all contain the same added barcode,
then we infer that they must have originated from the same original DNA molecule. Different barcode = different molecule of origin.
If the sequences would map to the same genomic region during sequence alignment, we would know that those sequences with the same
barcode originated from **a single DNA fragment from a single homologous chromosome from a single cell**. That's right, **built-in phase information**.
![simplified overview of linked reads](/static/linked_reads.svg)

### What do they look like?
Linked-read data is sequence data as you would expect it, encoded in a FASTQ file. The first processing step of
linked-read data is demultiplexing to split the raw Illumina-generated batch FASTQ file into samples (if multisample)
and identify/validate the linked-read barcode on every sequence. For 10X data, the barcode would _stay_ inline with
the sequence (to make it LongRanger compatible), but for other varieties (haplotagging, stLFR, etc.) you would also
remove the barcode from the sequence and preserve it using the `BX:Z` tag in the sequence headers. The demultiplexing process
is generally similar between non-10X linked-read technologies: a nucleotide barcode gets identified and moved from the sequence
to the read header under the `BX:Z` tag. The diagram below preserves the nucleotide barcode under the `OX:Z` tag and recodes
it under `BX:Z` using the haplotagging "ACBD" segment-aware format, however it would also be valid to just keep the nucleotide
barcode under `BX:Z`. Linked-read software is variable in its flexibility towards barcode formatting.
![haplotagging linked read data before and after demultiplexing](/static/lr_conversion.png)

#### Deconvolution
There are approaches to deconvolve linked read data (or "deconvolute", if you're indifferent to the burden of an extra syllable).
In this context, deconvolution is finding out which sequences are sharing barcodes by chance rather than because they originated
from the same DNA molecule. The likelihood of it happening is usually low, but it's not impossible. When the algorithm determines
that barcodes are being shared by unrelated molecules, the barcode typically gets suffixed with a hyphenated number (e.g.
`BX:Z:ATACG` becomes `BX:Z:ATACG-1`). As of this writing, there are only a few pieces of software that can deconvolve linked-read data and
do so with varying degrees of success and computational resource requirements. Similarly, linked-read software is variable in its
flexibility towards accepting the deconvolved-barcode format.

### Linked-read varieties
There are a handful of linked-read sample preparation methods ([read below](#linked-read-fastq-types)), but that's largely an implementation detail. All of those methods are
laboratory procedures to take genomic DNA and do the necessary modifications to fragment long DNA molecules, tag the resulting fragments with the same
DNA barcode, then add the necessary Illumina adapters. It's not unlike the different RAD flavors (e.g. EZrad, ddRAD, 2B-rad)-- they all give you RAD data in the end,
but vary in how you get there in terms of cost and bench time. We obviously subscribe to haplotagging :grin:.

### Failsafe
Unlike some other fancy well-touted sample preparation methods (_cough cough mate-pair_), linked-read data **is** whole genome
sequencing (WGS). What that means is that whether you use the linked-read information or not, the data will always be standard and 
viable WGS compatible with whatever you would use WGS for. It's WGS, but with a little extra info that goes a long way.

## Linked-read library performance
Because linked-read data has an extra dimension, there are some additional metrics that are useful to look at to evaluate how "good" 
the linked-read library turned out. As a baseline for comparison, think of RAD or WGS data and the kinds of metrics you might look at 
for library performance:
- coverage depth
- coverage breadth
- PCR/optical duplicates
- coverage depth per sample

Linked-read library performance also looks at that, but there's a few extra parameters that help us assess performance:

### Molecule Coverage
Since linked reads are tagger per DNA molecule, we are interested in understanding coverage on a per-molecule basis. By "molecule",
we are referring to the original DNA molecule from which barcoded sequences originated from.

#### reads per molecule
It's quite rare that all fragments of a single DNA molecule will get sequenced, so we are interested in getting an average number of
reads (sequences) per original DNA molecule. This is done by getting counts of all the reads with the same barcode, as unique barcodes
correspond with unique DNA molecules. It's difficult to say exactly what a good number of reads per molecule would be, but >2 would be
a good starting point (1 would be a singleton, see below). Having 3-6 would be decent, depending on your project goals.

#### percent molecule coverage
Because only a few fragments from a DNA molecule end up getting sequenced, it would be good to understand how much of a molecule is 
represented in seqences (breadth). It's impossible to get an accurate calculation for this metric because we have no way of knowing 
the size of the original DNA molecule, but what we can do is, after alignment, find the distance between the two furthest reads 
sharing a barcode and calculate how many bases between them have been sequenced. It's not perfect, but it's something. A low number 
would indicate large gaps between linked sequences (e.g. few sequences and/or a very large molecule), whereas a high number would 
indicate small gaps between linked sequences (e.g. many sequences and/or a very small molecule).
![molecule coverage example for 3 linked-reads sharing the same barcode](/static/molecule_coverage.svg)

### Singletons
Because linked reads need to be, well, _linked_, we need to know exactly how many of the sequences actually share barcodes. A barcode 
that only appears in one paired or unpaired read is a **singleton**, meaning it isn't actually linked to any other sequence. Despite 
having a linked-read barcode, the absence of other reads with the same barcode means the barcode information for that read (or read 
pair) is mostly useless, aka it's just a plain-regular short read and can be used as such. Having <60% would be considered decent, 
but you would want as few singletons as possible in a linked-read dataset. For perspective, if 90% of your reads were singletons, you 
can think of that as "90% of your data aren't linked reads". There isn't yet a consensus for what
to name the opposite, which is when a barcode is shared by more than one paired or unpaired read (a proper linked read), but our team 
calls them **nonsingletons**. It's admittedly clunky and we are open to suggestions.
![singletons and nonsingletons](/static/linked_singletons.svg)

## Linked Coverage
Because of the "linked" component of linked-read data, we have an additional kind of sequence coverage
to consider, which is **linked coverage**. Since linked-read barcodes preserve phase information, we can do an alternative kind of
coverage calculation where you pretend that all the sequences with a shared barcode are one big gapless sequence. If you use your 
imagination that way, the name "synthetic long reads" begins to make more sense. Mathematically, the linked coverage (breadth and 
depth) should always be higher than the coverage from just the alignments themselves, as gaps between linked reads are counted as 
well.
![linked depth](/static/linked_depth.svg)

In real data, this actually looks very cool. Below is a graph from the per-sample alignment report Harpy generates. It's a circos 
plot, which is a circular representation of round charts, and each wedge is a chromosome (start to finish), labelled by the "2R", 
"2L", etc. around the perimeter. The inner circle (grey histogram) is the **sequence alignment depth** (i.e. the standard depth 
calculation) and the outer circle (magenta histogram)
is the **linked depth** of those same data in those same 50kbp intervals. It has the most fascinating flower petal pattern, which
is due to the lower likelihood of reads linking at the edges of chromosomes.
![linked depth in real data](/static/linked_depth_example.png)

## Linked-Read FASTQ types
This isn't really the place to go into the nitty-gritty of the different linked-read chemistries
(i.e. I don't actually know the nuanced details), but it's worth describing the obvious differences
of the raw (FASTQ) data. Knowing these details might help you make sense of compatibilties/incompatibilities
for software, or how you can convert between styles.

### Haplotagging
Unlike the available chemistries, haplotagging is non-commercial (DIY!!). Haplotagging barcodes are combinatorial and
are made up of four 6bp segments. Two of these segments ("A" and "C") are the first 12bp of the I1 read and
the other two ("B" and "D") are the first 12bp of the I2 read, both of which are provided by Illumina for standard sequencing runs.
Because of this segment design, there are $96^4$ (~84 million) possible barcode combinations (~900,000 per sample). 
The barcodes are stored in the sequence header under the `BX:Z` SAM tag, recoded in their "`ACBD`" format.
- 4 barcode segments
  - `A` segment is the first 6bp of the I1 read
  - `C` segment is the next 6bp of the I1 read (7-12)
  - `B` segment is the first 6bp of the I2 read
  - `D` segment is the next 6bp of the I2 read (7-12)
- barcode stored as `BX:Z` tag in the read header in `ACBD` format
  - e.g. `@A003432423434:1:324 BX:Z:A45C01B84D21`

### TELLseq
One of the presently available commercial linked-read options. TELLseq data is very similar to 10X, except the
barcode is 18bp long and contained in the I1 read that Illumina provides with the standard R1 and R2 reads. The
barcode gets appended in the read header using a colon (`:`).
- barcode is the first 18bp of the I1 read 
- barcode is appended to sequence header
  - e.g. `@A00234534562:1:544:AATTATACCACAGCGGTA`
- advertised to have a capacity over 2 **billion** barcodes, but realistically use <24 **million**

### stLFR
Another of the presently available commercial linked-read options. stLFR data uses combinatorial barcodes
made up for three 10bp segments which are at the end of the R2 read. Demultiplexing these data results
in the barcode being moved to the sequence ID using a pound (`#`) sign between the sequence ID and barcode, with
the barcode recoded in the `1_2_3` format, where each segment is an integer.
- depending on the link sequence between segments, will be either the last 54bp or 42bp of the R2 read
  - 54 base barcode: 10+6+10+18+10
  - 42 base barcode: 10+6+10+6+10
- barcode appended to sequence header with `#` sign
  - e.g. `@A003432423434:1:324#12_432_1`
- advertised to have a capacity over 3.6 **billion**, with up to 50 **million** per sample (actual results may vary)

### 10X
The elder and thus least advanced of the bunch, and also a discontinued commercial product. 10X-style FASTQ files have the linked-read barcode
as the first 16bp of the forward (R1) read. For these data to be compatible with the 10X Longranger suite,
the barcode **must** stay in the read. Moving the first 16bp into the read header breaks Longranger compatibility.
- barcode is the first 16bp of the R1 read
- barcode stays in the sequence data for LongRanger compatibility
- limited to ~4.7 million barcodes

### Standard
We here at Harpy would like to _plead_ for linked-read practitioners to adopt a standardized representation of linked-read barcodes.
We don't really care what linked-read technology you use (although we obviously prefer haplotagging), but we do care about
the data formats being consistent. The classic joke in population genetics is that before you write a new piece of software,
you need to first write a new file format. We really don't want that here-- inconsistent linked-read data formats make it harder
for developers to write great software that is compatible with all linked-read data types when they need to account for every possible
location a barcode could be in, what a valid vs invalid barcode looks like for that technology, etc. So, we are humbly asking 
(nay, begging) for there to be a standard data format (FASTQ and SAM/BAM) that uses two SAM tags to encode barcodes:
- `BX:Z` tag to record the barcode, the format of which is irrelevant
  - it could be haplotagging `ACBD`, stLFR `1_2_3`, nucleotides, _whatever_
- `VX:i` tag to record if the corresponding barcode is valid
  - `0` is invalid
  - `1` is valid

This format makes it the responsibility of the early data processing software 
(as early as demultiplexing) to encode the data in this format for downstream
processing. Our hope is that this htslib-compliant format is adopted to reduce
the fragmentation in the software ecosystem and reduce the need for file format
conversions. It also creates a blueprint for a generic file encoding for any new
linked-read methods that are being developed or can be developed in the future.

#### FASTQ
``` stLFR style valid FASTQ header
@A00470:481:HNYFWDRX2:1:2101:29532:1063/1 BX:Z:45_11_361 VX:i:1
```
``` haplotagging style invalid FASTQ header
@A00470:481:HNYFWDRX2:1:2101:29532:1063/1 BX:Z:A78C00B14D96 VX:i:0
```
#### SAM/BAM
``` haplotagging style valid SAM record
A00470:481:HNYFWDRX2:1:2229:29912:29778 99      2R      1222    40      19S80M  =   1500     428     AGATGTGTATAAGAGACAGAGTTATGTCATTTTAAGCGGTCAAAATGGGTGAATTTCCGATTTCAAGTAATAGGCGAACTCAAGATACCTTCTACAGAT  FFFFFFFFFFFFFFFFF:FFFFF:FFFFF:FFFF:FFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:,FFFFFFF:FFFFF  NM:i:0  MD:Z:80 MC:Z:150M       AS:i:80      XS:i:80 RG:Z:sample1    MI:i:1  VX:i:1  BX:Z:A55C67B91D96
```

``` nucleotide (TELLseq/10X) style invalid SAM record
A00470:481:HNYFWDRX2:1:2229:29912:29778 99      2R      1222    40      19S80M  =   1500     428     AGATGTGTATAAGAGACAGAGTTATGTCATTTTAAGCGGTCAAAATGGGTGAATTTCCGATTTCAAGTAATAGGCGAACTCAAGATACCTTCTACAGAT  FFFFFFFFFFFFFFFFF:FFFFF:FFFFF:FFFF:FFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:,FFFFFFF:FFFFF  NM:i:0  MD:Z:80 MC:Z:150M       AS:i:80      XS:i:80 RG:Z:sample1    MI:i:1  VX:i:0  BX:Z:TACANNNCACAGAG
```