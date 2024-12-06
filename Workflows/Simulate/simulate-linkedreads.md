---
label: Linked Reads
description: Simulate linked reads from a genome
icon: dot
#visibility: hidden
order: 6
---

# :icon-flame: Simulate Linked Reads
Simulate linked reads from a genome

===  :icon-checklist: You will need
- two haplotypes of a reference genome in FASTA format: [!badge variant="success" text=".fasta"] [!badge variant="success" text=".fa"] [!badge variant="success" text=".fasta.gz"] [!badge variant="success" text=".fa.gz"] [!badge variant="secondary" text="case insensitive"]
    - can be created with [!badge corners="pill" text="simulate {snpindel,inversion,...}"](simulate-variants.md)
    - see the [tutorial](/blog/simulate_diploid.md)
- [!badge variant="ghost" text="optional"] a file of barcodes to tag linked reads with
==- :icon-question: LRSIM differences
The original [LRSIM](https://github.com/aquaskyline/LRSIM) is a lengthy Perl script that, like Harpy, outsources
to various other programs (`SURVIVOR`, `DWGSIM`, `samtools`, `msort`) and acts as a workflow through these programs. The Harpy
version of `LRSIM` keeps only the novel `LRSIM` code that creates linked reads from reads simulated by `DWGSIM`. The
rest of `LRSIM`'s components are reincorporated into the Snakemake workflow governing the [!badge corners="pill" text="simulate linkedreads"]
module, while removing the `SURVIVOR` part since [!badge corners="pill" text="simulate {snpindel,...}"](simulate-variants.md) are used for that purpose.
We generally wouldn't recommend using the Harpy version of LRSIM (`HaploSim.pl`) by itself unless you're confident you know what you're doing.
#### Notable differences
- dependencies are expected to be on the `PATH`, not hardcoded to the folder `LRSIM` is running from
- `-r`, `-u` parameters removed
- `-p` names the outputs with a prefix
- `-g` is used for the DWGSIM output fastq files, not input genomes
- `-a` is used for input fasta.fai files
- accepts any length of barcodes (mininmum 6 bases)
- **does not** add sequencing error to the barcode
- no `.status` file, all logging writes to `stderr`
- `SURVIVOR` variant simulation functionality removed entirely
- `DWGSIM`, samtools, msort, and extractReads functionality moved into Harpy workflow
- uses newer version of `DWGSIM`
===

You may want to benchmark haplotag data on different kinds of genomic variants. To
do that, you'll need *known* variants (like those created by  [!badge corners="pill" text="simulate {snpindel,...}"](simulate-variants.md)) and
linked-read sequences. This module will create (diploid) linked-read sequences from two genome haplotypes.
To accomplish this, Harpy uses a modified version of [LRSIM](https://github.com/aquaskyline/LRSIM),
and demultiplexes the resulting linked-reads into Haplotag-formatted reads. To simulate linked reads, use:

```bash usage
harpy simulate linkedreads OPTIONS... HAP1_GENOME HAP2_GENOME
```
```bash example
harpy simulate linkedreads -t 4 -n 2  -l 100 -p 50  data/genome.hap1.fasta data/genome.hap2.fasta
```

## :icon-terminal: Running Options
In addition to the [!badge variant="info" corners="pill" text="common runtime options"](/commonoptions.md), the  [!badge corners="pill" text="simulate linkedreads"] module is configured using these command-line arguments:

{.compact}
| argument            | short name |              default               | description                                                                                 |
| :------------------ | :--------: | :--------------------------------: | :------------------------------------------------------------------------------------------ |
| `HAP1_GENOME`       |            |                                    | [!badge variant="info" text="required"] Haplotype 1 of the diploid genome to simulate reads |
| `HAP2_GENOME`       |            |                                    | [!badge variant="info" text="required"] Haplotype 2 of the diploid genome to simulate reads |
| `--barcodes`        |    `-b`    | Meier et al. haplotagging barcodes | File of linked-read barcodes to add to reads                                                |
| `--distance-sd`     |    `-s`    |                `15`                | Standard deviation of read-pair distance                                                    |
| `--molecule-length` |    `-l`    |               `100`                | Mean molecule length (kbp)                                                                  |
| `--molecules-per`   |    `-m`    |                `10`                | Average number of molecules per partition                                                   |
| `--mutation-rate`   |    `-r`    |              `0.001`               | Random mutation rate for simulating reads (0 - 1.0)                                         |
| `--outer-distance`  |    `-d`    |               `350`                | Outer distance between paired-end reads (bp)                                                |
| `--patitions`       |    `-p`    |               `1500`               | Number (in thousands) of partitions/beads to generate                                       |
| `--read-pairs`      |    `-n`    |               `600`                | Number (in millions) of read pairs to simulate                                              |

## Mutation Rate
The read simulation is two-part: first `dwgsim` generates forward and reverse FASTQ files from the provided genome haplotypes
(`HAP1_GENOME` and `HAP2_GENOME`), then `LRSIM` takes over and creates linked-reads from that. The `--mutation-rate`
option controls random mutation rate `dwgsim` uses when creating FASTQ files from your provided genome haplotypes. This parameter
adds SNPs/variation in addition to the error rate assumed for the Illumina platform. If you don't want any more SNPs added to the
reads beyond sequencing error, set this value to `--mutation-rate 0`.

#### Simulating a single sample
If you intend to simulate a "single individual" (i.e. use this module once), then you might want no additonal SNPs beyond the variants
you may have already introduced into the genome and set `--mutation-rate 0`.

#### Simulating multiple samples
If you intend on simulating "multiple individuals" (i.e. use this module multiple times on the same genome haplotypes),
it may make sense to set this value larger than 0 so there is some "natural" variation between your simulated individuals.

## Partitions
**TL;DR**: 10X partitions ≈ haplotag beads

The option `--partitions` refers to the reaction "bubbles" in the original 10X linked-read chemistry. The 10X
protocol involved emulsion reactions where microscopic bubbles resulting from emulsion were each their own
reaction micro-environment. Each of these "partitions" (_aka_ bubbles, etc.) was to contain a unique linked-read
barcode that would be ligated onto the sample DNA, thus creating the linked read barcoding. In an ideal situation,
there would be a single molecule per emulsion partition, which was rarely the case because it's really
difficult to achieve that. In haplotag terms, think of partitions as being synonymous with tagmentation beads. In
both 10X and haplotag simulation cases, you're controlling for how many clashing barcodes there might be where
reads that aren't from the same molecule have the same linked-read barcode. This is why linked-read software (and
corresponding Harpy modules) have an option to set the [barcode threshold](../../haplotagdata.md/#barcode-thresholds). 

## Barcodes
Barcodes, if provided, must be given as nucleotide sequences, one per line. If not provided,
Harpy will generate the standard set of $96^4$ haplotagging barcodes (24bp) used by Meier et al.
**All barcodes must be the same length**. A barcode file should look like this:
``` input.barcodes.txt
ATATGTACTCATACCA
GGATGTACTCATTCCA
TCACGTACTCATACCA
etc...
```

### Haplotagging conversion
Harpy will convert (demultiplex) the simulated linked-reads into proper haplotagging ACBD format by matching the first `X`
number of nucleotides from the start of the forward read to the barcode list. The provided barcodes will all be assigned a
unique haplotagging ACBD barcode and moved to the `BX` tag in the sequence headers, e.g. `BX:Z:A01C47B32D91`. The original nucleotide basrcode 
will be removed from the sequence and stored in the `OX:Z` sequence header tag, e.g. `OX:Z:ATATGTACTCATACCA`.
The paired reverse read will also have these tags. The diagram below attempts to simplify this visually.
![inline linked-read barcode conversion into AxxCxxBxxDxx haplotag barcode format](/static/lr_conversion.png)

## Choosing parameters
LRSIM does internal calculations to determine the number of reads per molecule based on `--read-pairs`,
`--partitions`, and `--molecules-per`. Understanding how these parameters affect the resulting sequences
will help inform your decisions for those parameters:

$$
\text{Reads Per Molecule} = 0.499 + \frac{N \times 1,000,000}{\left(\frac{P \times 1,000}{H}\right) \times M \times H}
$$
$$\text{where:}\\\text{N = number of reads to simulate (in millions)}\\\text{H = number of haplotypes (fixed at 2)}\\\text{P = number of partitions (in thousands)}\\\text{M = molecules per partition}$$

### Parameter calculator
Use the calculator provided below to help you make informed decisions for these parameters:
[!embed](https://app.calconic.com/api/embed/calculator/662146310482ea001e7acea2)

## :icon-git-pull-request: Simulate Linkedreads Workflow

```mermaid
graph LR
    subgraph Inputs
        direction BT
        A[genome haplotype 1]:::clean
        B[genome haplotype 2]:::clean
    end
    Inputs-->D([dwgsim]):::clean
    D-->L([LRSIM]):::clean
    L-->H([convert to haplotag]):::clean
    style Inputs fill:#f0f0f0,stroke:#e8e8e8,stroke-width:2px
    classDef clean fill:#f5f6f9,stroke:#b7c9ef,stroke-width:2px
```
