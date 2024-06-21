---
label: Align
description: Align haplotagged sequences
---
# :icon-quote: Align Sequences to a Genome

After your sequences (in FASTQ format) have been checked for quality, you
will need to align them to a reference genome before you can call variants.
Harpy offers several aligners for this purpose:

{.compact}
| aligner | linked-read aware | speed | link |
| :--- | :---: | :---:| :---: |
| [BWA](bwa.md) | no ❌ | fast ⚡ | [repo](https://github.com/lh3/bwa), [paper](http://arxiv.org/abs/1303.3997) |
| [EMA](ema.md) | yes ✅ | slow 🐢 |[repo](https://github.com/arshajii/ema), [paper](https://www.biorxiv.org/content/early/2017/11/16/220236) |
| [Minimap2](minimap.md) | no ❌ | fast ⚡ | [repo](https://github.com/lh3/minimap2) [paper](https://doi.org/10.1093/bioinformatics/btab705) |

Despite the fact that EMA is the only barcode-aware aligner offered, when using BWA or Minimap2, Harpy retains the barcode information from the sequence headers and will
assign molecule identifiers (`MI:i` SAM tags) based on these barcodes and the [molecule distance threshold](../../haplotagdata.md/#barcode-thresholds).