---
label: Downsample
description: Downsample data by barcode
icon: fold-down
order: 10
---

# :icon-fold-down: Downsample data by barcode

===  :icon-checklist: You will need one of either
- one alignment file [!badge variant="success" text=".bam"] [!badge variant="success" text=".sam"] [!badge variant="secondary" text="case insensitive"]
- one set of paired-end reads in FASTQ format [!badge variant="success" text=".fq"] [!badge variant="success" text=".fastq"] [!badge variant="secondary" text="gzip recommended"] [!badge variant="secondary" text="case insensitive"]
===

While downsampling (subsampling) FASTQ and BAM files is relatively simple with tools such as `awk`, `samtools`, `seqtk`, `seqkit`, etc.,
[!badge corners="pill" text="downsample"] allows you to downsample a BAM file (or paired-end FASTQ) _by barcodes_. That means you can
keep all the reads associated with `d` number of barcodes. The `--invalid` proportion will determine what proportion of invalid barcodes appear in the barcode
pool that gets subsampled, where `0` is none, `1` is all invalid barcodes, and a number in between is that proportion, e.g. `0.5` is half.
Bear in mind that the barcode pool still gets subsampled, so the `--invalid` proportion doesn't necessarily reflect how many end up getting
sampled, rather what proportion will be considered for sampling. 

!!! Barcode tag
Barcodes must be in the `BX:Z` SAM tag for both BAM and FASTQ inputs. See [Section 1 of the SAM Spec here](https://samtools.github.io/hts-specs/SAMtags.pdf).
!!!

```bash usage
harpy downsample OPTIONS... INPUT(S)...
```

```bash example
# BAM file
harpy downsample -d 1000 -i 0.3 -p sample1.sub1000 sample1.bam

# FASTQ file
harpy downsample -d 1000 -i 0 -p sample1.sub1000 sample1.F.fq.gz sample1.R.fq.gz
```

## :icon-terminal: Running Options
In addition to the [!badge variant="info" corners="pill" text="common runtime options"](/commonoptions.md), the [!badge corners="pill" text="downsample"]
module is configured using the command-line arguments below.

{.compact}
| argument        | short name |    default    | description                                                                                                                       |
| :-------------- | :--------: | :-----------: | :-------------------------------------------------------------------------------------------------------------------------------- |
| `INPUT(S)`      |            |               | [!badge variant="info" text="required"] One BAM file or both read files from a paired-end FASTQ pair                              |
| `--downsample`  |    `-d`    |               | [!badge variant="info" text="required"] Number of barcodes to downsample to                                                       |
| `--invalid`     |    `-i`    |      `1`      | Proportion of barcodes to sample                                                                                                  |
| `--prefix`      |    `-p`    | `downsampled` | Prefix for output files                                                                                                           |
| `--random-seed` |            |               | Random seed for sampling [!badge variant="secondary" text="optional"]                                                             |

----
## :icon-git-pull-request: Downsample Workflow
```mermaid
graph LR
    subgraph fastq
        R1([read 1]):::clean---R2([read 2]):::clean
    end
    subgraph bam
        bamfile([bam]):::clean
    end
    fastq-->|bam conversion|bam
    bam-->sub([extract and\n subsample barcodes]):::clean
    sub-->exreads([extract reads]):::clean
    bam-->exreads
    fastq-->exreads
    style fastq fill:#f0f0f0,stroke:#e8e8e8,stroke-width:2px
    style bam fill:#f0f0f0,stroke:#e8e8e8,stroke-width:2px
    classDef clean fill:#f5f6f9,stroke:#b7c9ef,stroke-width:2px
```