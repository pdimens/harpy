---
label: Align
description: Align haplotagged sequences
---
# :icon-quote: Align Sequences to a Genome

After your sequences (in FASTQ format) have been checked for quality, you
will need to align them to a reference genome before you can call variants.
Harpy offers several aligners for this purpose:

{.compact}
| aligner | linked-read aware | speed | repository | publication |
| :--- | :---: | :---:| :---: | :---:|
| [BWA](bwa.md) | no ❌ | fast ⚡ | [github](https://github.com/lh3/bwa) | [paper](http://arxiv.org/abs/1303.3997) |
| [EMA](ema.md) | yes ✅ | slow 🐢 |[github](https://github.com/arshajii/ema) | [preprint](https://www.biorxiv.org/content/early/2017/11/16/220236) |
| [strobealign](strobe.md) | no ❌ | super fast ⚡⚡ | [github](https://github.com/ksahlin/strobealign) | [paper](https://doi.org/10.1186/s13059-022-02831-7) |

Despite the fact that EMA is the only barcode-aware aligner offered, when using BWA or strobealign, Harpy retains the barcode information from the sequence headers and will
assign molecule identifiers (`MI:i` SAM tags) based on these barcodes and the [molecule distance threshold](../../haplotagdata.md/#barcode-thresholds).

## [!badge text="New"] Works with regular WGS data
Starting with Harpy `v2.x`, the `--ignore-bx` option lets you skip the workflow
routines that do things specific to linked reads, meaning you can comfortably use
[!badge corners="pill" text="harpy align bwa"](bwa.md) and [!badge corners="pill" text="harpy align strobe"](strobe.md) to align your WGS sequence data. 

!!!warning RADseq data
RADseq data will probably work fine too, however you may need to post-process the
BAM files to unset the duplicate flag, as marking duplicates in RADseq (without UMIs) [may cause issues](https://www.researchgate.net/post/How_to_exclude_PCR_duplicates_in_ddRAD) with SNP calling:
```bash
samtools view -b -h --remove-flags 1024 -o output.bam input.bam
```
!!!