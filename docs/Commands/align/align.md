---
label: align
description: Choosing an aligner
---
# :icon-quote: Align Sequences to a Genome

After your sequences (in FASTQ format) have been checked for quality, you
will need to align them to a reference genome before you can call variants.
Harpy offers several aligners for this purpose:

{.compact .clean .whitespace-nowrap}
| command | aligner     | best for          |                                       repository |                      publication                       |
| :------ | :---------- | :---------------- | -----------------------------------------------: | :----------------------------------------------------: |
| bwa     | minibwa     | general use       |         [github](https://github.com/lh3/minibwa) |       [preprint](https://github.com/lh3/minibwa)       |
| minimap | minimap2    | long reads        |        [github](https://github.com/lh3/minimap2) | [paper](https://doi.org/10.1093/bioinformatics/bty191) |
| strobe  | strobealign | speed             | [github](https://github.com/ksahlin/strobealign) |  [paper](https://doi.org/10.1186/s13059-022-02831-7)   |
| arachne | arachne     | linked-read aware |    [github](https://github.com/blinkseq/arachne) |                                                        |

Neither minibwa, minimap2, nor strobealign are linked-read aware aligners, but Harpy transfers the barcode information from the sequence headers into the alignments.
Use `arachne` (experimental) for using linked-read information to resolve alignments over repeat regions. 

==- Linked-read alignment
Arachne is a fork and successor of Lariat, the aligner written by 10X Genomics for the original linked-read technology.
It defaults to using an improved version of the EM algorithm from [EMA](https://www.biorxiv.org/content/10.1101/220236v1),
but can use the RFA algorithm from Lariat with the `--rfa` flag. Incorporating linked-read information during mapping
improves accuracy over barcode-naive alignment, as shown here:

![EMA Publication figure 3](/static/EMA.fig3.png)
===

## Non linked-read WGS data
Starting with Harpy `v2.x`, you can skip the workflow
routines that do things specific to linked reads, meaning you can comfortably use
[!badge corners="pill" text="harpy align bwa"](standard.md) and [!badge corners="pill" text="harpy align strobe"](standard.md) to align your WGS sequence data. 
- version `2.0-2.7` : `--ignore-bx`
- version `>2.7` : `--lr-type none`
- version `>=3.0`: autodetected or forced with `--unlinked`

!!!warning RADseq data
RADseq data will probably work fine too, however you may need to post-process the
BAM files to unset the duplicate flag, as marking duplicates in RADseq (without UMIs) [may cause issues](https://www.researchgate.net/post/How_to_exclude_PCR_duplicates_in_ddRAD) with SNP calling:
```bash
samtools view -b -h --remove-flags 1024 -o output.bam input.bam
```
!!!
