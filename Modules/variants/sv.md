---
label: SV
description: Call structural variants on alignments generated from haplotagged sequences with Harpy
icon: project-roadmap
order: 1
---

# :icon-sliders: Call Structural Variants
(like indels, insertions, duplications)

===  :icon-checklist: You will need
- at least 4 cores/threads available
- a genome assembly in FASTA format
- sequence alignments, in `.bam` format
- sample grouping file ([see below](#pooled-sample-variant-calling))
===

After reads have been aligned, _e.g._ with `harpy align`, you can use those alignment files
(`.bam`) to call variants in your data. Harpy can call structural variants using [LEVIATHAN](#leviathan-workflow) or [NAIBR](#naibr-workflow),
which identify inversions, duplications, and deletions. You can call variants with Harpy using the
`variants sv` module:

```bash usage
harpy variants sv OPTIONS... 
```

```bash examples
# call structural variants with naibr
harpy variants sv --threads 20 --genome genome.fasta --directory Align/bwa

# call structural variants with LEVIATHAN
harpy variants sv --threads 20 --genome genome.fasta --directory Align/bwa --method leviathan
```

## :icon-terminal: Running Options
In addition to the [common runtime options](../../commonoptions.md), the `harpy variants sv` module is configured using these command-line arguments:

| argument         | short name | type                            | default | required | description                                        |
|:-----------------|:----------:|:--------------------------------|:-------:|:--------:|:---------------------------------------------------|
| `--genome`       |    `-g`    | file path                       |         | **yes**  | Genome assembly for variant calling                |
| `--directory`    |    `-d`    | folder path                     |         | **yes**  | Directory with sequence alignments                 |
| `--populations`  |    `-p`    | file path                       |         |    no    | Tab-delimited file of sample\<*tab*\>group           |
| `--ploidy`       |    `-x`    | integer                         |    2    |    no    | Ploidy of samples                                  |
| `--method`       |    `-l`    | choice [`naibr`, `leviathan`] | mpileup |    no    | Which variant caller to use                          |
| `--extra-params` |    `-x`    | string                          |         |    no    | Additional naibr/leviathan arguments, in quotes    |

### The --populations option
#### Single-sample variant calling
When **not** using a population grouping file via `--populations`, variants will be called per-sample. 
Due to the nature of Structural Variant (SV) VCF files, there isn't an entirely fool-proof way 
of combining the variants of all the samples into a single VCF file, therefore the output will be a VCF for every sample.

#### Pooled-sample variant calling
With the inclusion of a population grouping file via `--populations`, Harpy will merge the bam files of all samples within a 
population and call SV's on these alignment pools. Preliminary work shows that this way identifies more variants and with fewer false 
positives. **However**, individual-level information gets lost using this approach, so you will only be able to assess 
group-level variants, if that's what your primary interest is. 

|||:icon-file: sample grouping file
This file is optional and only useful if you want variant calling to happen on a per-population level.
- takes the format of sample\<*tab*\>group
- the groups can be numbers or text (_i.e._ meaningful population names)
- create with `harpy extra -p <samplefolder>` or manually
- if created with `harpy extra -p`, all the samples will be assigned to group `pop1`, so make sure to edit the second column to reflect your data correctly.

``` example file for --populations
sample1 pop1
sample2 pop1
sample3 pop2
sample4 pop1
sample5 pop3
```

!!!warning known quirk
There's an unusual error on the Snakemake side of things that happens when the name of a sample and population are identical.
It has been unclear how to resolve this issue, so to protect yourself, it's best to make sure the population names are different
from the sample names. A simple fix would be to use underscores (`_`) to differentiate the population name.
!!!
|||
----

### :icon-git-pull-request: NAIBR workflow
+++ :icon-git-merge: details
[Naibr](https://github.com/raphael-group/NAIBR) is a variant caller that uses linked read barcode information 
to call structural variants (indels, inversions, etc.) exclusively, meaning it does not call SNPs. The authors of Naibr have not been updating or improving it, so Harpy uses
[an active fork](https://github.com/pontushojer/NAIBR) of it that is available on [Bioconda](https://anaconda.org/bioconda/naibr-plus) under the name `naibr-plus`. This fork includes improved accuracy as well as quality-of-life updates.

```mermaid
graph LR
    subgraph Population calling
    popsplit([merge by population])
    end
    subgraph Individual calling
    bams([individual alignments])
    end
    popsplit-->A
    bams-->A
    A([index alignments]) --> B([naibr])
    Z([create config file]) --> B
    popsplit --> Z
    bams --> Z
    B-->C([generate reports])
```
+++ :icon-file-directory: naibr output
The `harpy variants --method naibr` module creates a `Variants/naibr` (or `naibr-pop`) directory with the folder structure below. `sample1` and `sample2` are generic sample names for demonstration purposes.

```
Variants/naibr/
├── sample1.bedpe
├── sample2.bedpe
├── configs
│   ├── sample1.config
│   └── sample2.config
├── filtered
│   ├── sample1.fail.bedpe
│   └── sample2.fail.bedpe
├── IGV
│   ├── sample1.reformat.bedpe
│   └── sample2.reformat.bedpe
├── logs
│   ├── harpy.variants.log
│   ├── sample1.log
│   └── sample2.log
├── reports
│   ├── sample1.naibr.html
│   └── sample2.naibr.html
└── vcf
    ├── sample1.vcf
    └── sample2.vcf
```

| item          | description                                                      |
|:--------------|:-----------------------------------------------------------------|
| `*.bedpe`     | structural variants identified by NAIBR                          |
| `configs/`    | the configuration files harpy generated for each sample          |
| `filtered/`   | the variants that failed NAIBR's internal filters                |
| `IGV/`        | same as the output .bedpe` files but in IGV format               |
| `logs/harpy.variants.log` | relevant runtime parameters for the variants module  |
| `logs/*.log`  | what NAIBR writes to `stderr` during operation                   |
| `reports/`    | summary reports with interactive plots of detected SV            |
| `vcf/`        | the resulting variants, but in `.VCF` format                     |

+++ :icon-code-square: naibr parameters
By default, Harpy runs `naibr` with these parameters (excluding inputs and outputs):
```python
min_mapq = 30
d        = 10000
min_sv   = 1000
k        = 3
```

Below is a list of all `NAIBR` runtime options, excluding those Harpy already uses or those made redundant by Harpy's implementation of NAIBR.
These are taken directly from the [NAIBR documentation](https://github.com/pontushojer/NAIBR#running-naibr). If adding these arguments, do so like:
`-x "min_sv 1000 d 50000"`
``` NAIBR arguments
 -d: The maximum distance in basepairs between reads in a linked-read (default: 10000)
 -blacklist: BED-file with regions to be excluded from analysis
 -candidates: BEDPE-file with novel adjacencies to be scored by NAIBR. This will override automatic detection of candidate novel adjacencies
 -min_sv: Minimum size of a structural variant to be detected (default: lmax, i.e. the 95th percentile of the paired-end read insert size distribution)
 -k: minimum number of barcode overlaps supporting a candidate NA (default = 3)
```
+++

### :icon-git-pull-request: LEVIATHAN workflow
+++ :icon-git-merge: details
[Leviathan](https://github.com/morispi/LEVIATHAN) is an alternative variant caller that uses linked read barcode information 
to call structural variants (indels, inversions, etc.) exclusively, meaning it does not call SNPs. Harpy first uses [LRez](https://github.com/morispi/LRez) to index the barcodes 
in the alignments, then it calls variants using Leviathan.

!!!warning EMA-mapped reads
Leviathan relies on split-read information in the sequence alignments to call variants. The EMA aligner
does not report split read alignments, instead it reports secondary alignments. It is recommended to use
BWA-generated alignments if intending to call variants with leviathan. 
!!!

```mermaid
graph LR
    subgraph Population calling
    popsplit([merge by population])
    end
    subgraph Individual calling
    bams([individual alignments])
    end
    popsplit-->A
    bams-->A
    A([index barcodes]) --> B([leviathan])
    B-->C([convert to BCF])
    C-->E([generate reports])
```
+++ :icon-file-directory: leviathan output
The `harpy variants --method leviathan` module creates a `Variants/leviathan` (or `leviathan-pop`) directory with the folder structure below. `sample1` and `sample2` are generic sample names for demonstration purposes.

```
Variants/leviathan/
├── sample1.bcf
├── sample2.bcf
├── logs
│   ├── harpy.variants.log
│   ├── sample1.leviathan.log
│   ├── sample1.candidates
│   ├── sample2.leviathan.log
│   └── sample2.candidates
├── reports
│   ├── sample1.SV.html
│   └── sample2.SV.html
└── stats
    ├── sample1.sv.stats
    └── sample2.sv.stats
```

| item                   | description                                              |
|:-----------------------|:---------------------------------------------------------|
| `*.bcf`                | structural variants identified by LEVIATHAN              |
| `logs/*.leviathan.log` | what LEVIATHAN writes to `stderr` during operation       |
| `logs/*candidates`     | candidate structural variants LEVIATHAN identified       |
| `reports/`             | summary reports with interactive plots of detected SV    |
| `logs/harpy.variants.log` | relevant runtime parameters for the variants module  |
| `stats/`               | results of `bcftools stats` on the vcf LEVIATHAN creates |

+++ :icon-code-square: leviathan parameters
By default, Harpy runs `leviathan` with default parameters (shown below), only modifying inputs and outputs at the command line.

Below is a list of all `leviathan` command line options, excluding those Harpy already uses or those made redundant by Harpy's implementation of LEVIATHAN.
These are taken directly from the [LEVIATHAN documentation](https://github.com/morispi/LEVIATHAN).
``` LEVIATHAN arguments
  -r, --regionSize:         Size of the regions on the reference genome to consider (default: 1000)
  -v, --minVariantSize:     Minimum size of the SVs to detect (default: same as regionSize)
  -n, --maxLinks:           Remove from candidates list all candidates which have a region involved in that much candidates (default: 1000)
  -M, --mediumSize:         Minimum size of medium variants (default: 2000)
  -L, --largeSize:          Minimum size of large variants (default: 10000)
  -s, --smallRate:          Percentile to chose as a threshold in the distribution of the number of shared barcodes for small variants (default: 99)
  -m, --mediumRate:         Percentile to chose as a threshold in the distribution of the number of shared barcodes for medium variants (default: 99)
  -l, --largeRate:          Percentile to chose as a threshold in the distribution of the number of shared barcodes for large variants (default: 99)
  -d, --duplicates:         Consider SV as duplicates if they have the same type and if their breakpoints are within this distance (default: 10)
  -s, --skipTranslocations: Skip SVs that are translocations (default: false)
  -p, --poolSize:           Size of the thread pool (default: 100000)
  -B, --nbBins:             Number of iterations to perform through the barcode index (default: 10)
  -c, --minBarcodes:        Always remove candidates that share less than this number of barcodes (default: 1)
```
+++