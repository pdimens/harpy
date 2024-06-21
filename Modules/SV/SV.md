---
label: SV
description: Find structural variants
---

# :icon-sliders: Find structural variants
The [!badge corners="pill" text="snp"](../snp.md) module identifies single nucleotide
polymorphisms (SNP) and small indels, but you may want to (and should!) leverage the linked-read
data to identify larger structural variants (SV) like large deletions, duplications, and
inversions. Harpy provides two linked-read variant callers to do exactly that:

{.compact}
| aligner | recommended | variants it detects | caveats | link |
| :--- | :---: | :--- | :--- |
| [NAIBR](naibr.md) |  ✅ | inversion, duplication, deletion | requires phased alignments | [original repo](https://github.com/raphael-group/NAIBR), [fork repo](https://github.com/pontushojer/NAIBR), [paper](https://doi.org/10.1093/bioinformatics/btx712) |
| [LEVIATHAN](leviathan.md) | | inversion, duplication, deletion, breakend | requires supplementary alignments |[repo](https://github.com/morispi/LEVIATHAN), [paper](https://doi.org/10.1101/2021.03.25.437002) |

## Caveats
### NAIBR
While our testing shows that NAIBR tends to find known inversions that LEVIATHAN misses, the program requires haplotype 
**phased** bam files as input. That means the alignments have a `PS` or `HP` tag that indicate
which haplotype the read/alignment belongs to. If your alignments don't have phasing tags (none of the current aligners in Harpy do this),
then you will need to use the [!badge corners="pill" text="phase"](../phase.md) module to phase your SNPs into haplotypes so
the [!badge corners="pill" text="sv naibr"](naibr.md) module can use that to phase your input alignments and proceed as planned.

### LEVIATHAN
LEVIATHAN relies on split-read information in the sequence alignments to call variants.
The EMA aligner does not report split read alignments, instead it reports secondary alignments.
It is recommended to use BWA- or Minimap2-generated alignments if intending to call variants with [!badge corners="pill" text="sv leviathan"](leviathan.md).