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
- at least one haplotypes of a reference genome in FASTA format: [!badge variant="success" text=".fasta"] [!badge variant="success" text=".fa"] [!badge variant="success" text=".fasta.gz"] [!badge variant="success" text=".fa.gz"] [!badge variant="secondary" text="case insensitive"]
    - can be created with [!badge corners="pill" text="simulate {snpindel,inversion,...}"](simulate-variants.md)
- to read [the Mimick documentation](https://pdimens.github.io/mimick/#/usage)
- [!badge variant="ghost" text="optional"] a file of barcodes to tag linked reads with
===

You may want to benchmark haplotag data on different kinds of genomic variants. To
do that, you'll need *known* variants (like those created by  [!badge corners="pill" text="simulate {snpindel,...}"](simulate-variants.md)) and
linked-read sequences. In Harpy `v1.x` this was done using a modified version of
[LRSIM](https://github.com/aquaskyline/LRSIM), however, Harpy `v2.x` now uses the purpose-built software [Mimick](https://github.com/pdimens/mimick)
(originally XENIA from the VISOR project). Mimick does exactly what you would need it to do, so
to keep the familiarity of [!badge corners="pill" text="simulate linkedreads"] in Harpy, we just expose a very thinly
veiled wrapper for Mimick. The only additions here are that Harpy automatically installs Mimick and can dispatch the job to
an HPC using `--hpc` like other workflows and ensures proper pairing and haplotype concatenation, otherwise it really just runs Mimick exactly as you would. All of Mimick's 
command line arguments are exposed to Harpy, except `--mutations`, `-indels`, and `-extindels`, which are set to `0`
to make sure you are only simulating linked-reads exactly.

Rather than having to maintain two copies of the same documentation, please
head over to [the Mimick documentation](https://pdimens.github.io/mimick/#/usage). 

```bash usage
harpy simulate linkedreads OPTIONS... BARCODES FASTA...
```
```bash example
harpy simulate linkedreads -t 4 18,96 data/genome.hap1.fasta data/genome.hap2.fasta
```

## :icon-terminal: Running Options
{.compact}
| long name                  | default         | type    | description                                                                                                                                                                                 |
|:---------------------------|:----------------|:--------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `BARCODES`                 |                 |         | [!badge variant="secondary" text="required"] input barcode file or length,count                                                                                                             |
| `FASTA`                    |                 |         | [!badge variant="secondary" text="required"] input fasta file(s)                                                                                                                            |
| `--output-prefix` `-o`     | `simulated/SIM` |         | output file prefix                                                                                                                                                                          |
| `--output-type` `-O`       | `standard`      |         | output format of FASTQ files                                                                                                                                                                |
| `--quiet` `-q`             | `0`             |         | `0` all output, `1` single progress bar, `2` no output                                                                                                                                      |
| `--regions` `-r`           |                 |         | one or more regions to simulate, in BED format                                                                                                                                              |
| `--threads` `-t`           | `2`             |         | number of threads to use for simulation                                                                                                                                                     |
| `--coverage`               | `30`            |         | mean coverage target for simulated data                                                                                                                                                     |
| `--distance`               | `500`           | integer | outer distance between the two ends in bp                                                                                                                                                   |
| `--error`                  | `0.02`          | `0`-`1` | base error rate, fixed at this number                                                                                                                                                       |
| `--extindels`              | `0`             | `0`-`1` | [!badge variant="secondary" text="hidden"] indel extension rate                                                                                                                             |
| `--indels`                 | `0`             | `0`-`1` | [!badge variant="secondary" text="hidden"] indel creation rate                                                                                                                              |
| `--length`                 | `150`           | `30`    | length of reads in bp, must be >`30`                                                                                                                                                        |
| `--mutation`               | `0`             | `0`-`1` | [!badge variant="secondary" text="hidden"] mutation rate                                                                                                                                    |
| `--stdev`                  | `50`            |         | standard deviation of `--distance`                                                                                                                                                          |
| `--lr-type` `-l`           | `haplotagging`  |         | type of linked-read experiment                                                                                                                                                              |
| `--molecule-coverage` `-c` | `0.2`           |         | mean percent coverage per molecule if <1, else mean number of reads per molecule                                                                                                            |
| `--molecule-length` `-m`   | `80000`         |         | mean length of molecules in bp                                                                                                                                                              |
| `--molecule-number` `-n`   | `3`             |         | mean number of unrelated molecules per barcode, where a negative number (e.g. `-2`) will use a fixed number of unrelated molecules and a positive one will draw from a Poisson distribution |

### Barcodes
#### Randomly generate
Mimick lets you put in length and count parameters, which it will use to randomly generate barcodes.
The format is `length,count` (no spaces), where `length` is the base-pair length the barcodes should be and
`count` is how many barcodes it should generate of length `length`. For example, if you specify `16,4000000`,
Mimick will generate 4 million unique 16bp barcodes, effectively mimicking 10X barcodes.
Mimick will write a file containing the barcodes it generated. In practice, this would look something like:
```randomly generated barcodes
harpy simulate linkedreads --lr-type 10x 16,4000000 hap1.fasta hap2.fasta
#                                     16bp^ ^4 million barcodes
```

#### Specific barcodes
Alternatively, if you have a set of barcodes you absolutely want to use, just put the filename as the first positional argument.
In practice, this would look something like:
```barcodes as a file
harpy simulate linkedreads --lr-type 10x barcodes.txt hap1.fasta hap2.fasta
#                                        ^file of barcodes
```

### FASTA file(s)
You will need at least 1 fasta file as input, which goes at the very end of the command. It's assumed
that each fasta file is a different haplotype **of the same genome**.
```fasta inputs
harpy simulate linkedreads --lr-type 10x 16,4000000 hap1.fasta hap2.fasta
#                                         ^barcodes  ^fasta 1   ^fasta 2
```
