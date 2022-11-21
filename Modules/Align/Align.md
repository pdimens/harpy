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

!!! Only aligned sequences
The aligners in this module only output aligned sequences, meaning unmapped reads will not appear in the results. Keep that in mind when you read the alignment summary report and see "100%" aligned reads. That just means there are no unmapped reads in the BAM file, which is the result of aligner output configuration and not biological processes.
!!!

Despite the fact that EMA is the only barcode-aware aligner offered, when using BWA or strobealign, Harpy retains the barcode information from the sequence headers and will
assign molecule identifiers (`MI:i` SAM tags) based on these barcodes and the [molecule distance threshold](../../haplotagdata.md/#barcode-thresholds).
