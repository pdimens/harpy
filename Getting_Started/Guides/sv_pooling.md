---
author: 
  - name: Pavel Dimens
    avatar: https://cals.cornell.edu/sites/default/files/styles/faculty/public/2024-09/afs-headshot-high-res-2cropped_0.jpg
date: 2024-08-01
category: guides
description: Why pool samples for SV calling and when to do it
icon: sync
label: SV calling pooled samples
image: https://visualpharm.com/assets/214/Merge%20Files-595b40b75ba036ed117d8636.svg
---

# :icon-sync: Pooling samples for SV calling
## TL;DR
- Samples can be concatenated together into a single "pool" and have SV called on the pool
- Pooling samples requires barcode deconvolution

One of the cool benefits of linked-read data is the fact that you
can call structural variants with it. Depending on the depth of 
your data, you may want (or need) to pool samples together. This 
page walks you through why you may want to do that and the logic of
doing so.

## Sample depth
### Depth, explained
In bioinformatics, the terms "coverage" and "depth" and often used
interchangeably, which is incorrect and leads to confusion. _Coverage_
refers to the proportion of a genome that is sequenced, and _depth_
refers to the number of sequences present at a given genomic position.
Think of coverage as horizontal and depth as vertical. If depth is given
as "15X", that means every sequenced genomic position has about 15 sequences
per position (a.k.a. 15 replicates per locus). For clarity, when we say "depth",
we are referring to this "vertical" description of replicates per locus.
![The difference between depth and coverage. The locus on the left would be considered 5X.](/static/depth_coverage.png)
### Depth, in context
Historically, one would have wanted to sequence fewer individuals at higher depth
to get confident genotype calls, rather than sequence more individuals at lower depth.
Recent advances in bioinformatics have enabled low-coverage whole genome sequencing
a.k.a. _lcWGS_ to be a viable option for studies of outbred non-model systems. More
and more studies are emerging that have adopted lcWGS, which is pretty cool. However,
calling structural variants at low depth is challenging, especially with short reads.

## The problem
It's recommended to have at least 10X-12X depth to get decent structural variant calls
(definitely read that in a paper that I would like to link here, but I can't seem to find
it). If your data already has a minimum of 10X for each individual, great! Feel free to use
variant callers like `naibr` and `leviathan` to identify structural variants. However, if
you opted to sequence more individuals at lower depth (lcWGS is often between 0.5-5X),
then calling structural variants in individuals may be a challenge.

## The solution
One way to get your low-coverage (low depth) data and still call structural variants is to pool
samples together, which would effectively boost the depth. By doing this, you will
no longer be able to make per-individual assessments, which can be fine depending on
the nature of your study.

## Pooling considerations
If pooling samples, you must pool them sensibly and with a biological context to do so.
In other words, you don't just pool random samples together to inflate depth. Since
haplotag data is just whole genome sequence data plus a little extra information, you should
use the SNPs of your data to first identify genetic clusters as per standard population
genetic practices. Once population structure/stratification has been identified, you can use
that as a basis to pool together samples from these groups. As an example, you can use
[pcangsd](https://github.com/Rosemeis/pcangsd) to perform a PCA on low-depth SNP data and
get results like these:
![PCA of Alosa sapidissima (SNPs from low-depth haplotag dataset)](/static/pca.png)

Given these results, a sensible pooling strategy may be:

{.compact}
| Pool | Sample Locations           |
| :---:| :------------------------- |
| 1    | Miramichi                  |
| 2    | Annapolis                  |
| 3    | St_Johns                   |
| 4    | Santee_Cooper              |
| 5    | Roanoke, Chowan-Blackwater |
| 6    | Potomac                    |
| 7    | Delaware                   |
| 8    | Hudson                     |
| 9    | Connecticut                |
| 10   | Merrimack, Kennebec        |

!!!danger Do not concatenate manually!
It would seem natural to use `samtools cat` to combine alignment files
quickly and easily, but **do not do this**. The [!badge corners="pill" text="harpy sv"](/Workflows/SV/SV.md)
workflows will intelligently concatenate files and will make sure
every individual will have unique `MI` values that are not shared with any
other individual in the pool. If you need to concatenate linked-read alignment files outside
of a workflow, use [concatenate_bam.py](/Getting_Started/Resources/utilities.md#concatenate_bam) shipped with Harpy
 instead of `samtools cat` or other similar tools. 
==- technical explanation
The reason is, samples aligned
using Harpy have their linked-read barcodes deconvolved and supplemented with
a unique molecule ID, given as an `MI:i` SAM tag in the alignment records. The
molecule ID's start at `1` for each sample and increment for every identified
unique molecule. What this means is that doing a typical concatenation via
`samtools cat` **will not resolve conflicting molecule IDs**.

Think of it this way,
every sample has MI's from `1` to `N`, and as soon as you do a basic concatenation,
the MI 1..N from sample2 will immediately conflict with 1..N from sample1. Why?
Well, it's impossible that `MI:i:1` from sample1 and `MI:i:1` from sample2 came from
the same molecule. A basic concatenation will flood the resulting file with clashing `MI`
tags, which will make it impossible for a linked-read aware SV caller to make sense
of the data and do its job well.
===
!!!
