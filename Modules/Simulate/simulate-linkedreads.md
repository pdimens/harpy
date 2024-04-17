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
- two haplotypes of a reference genome in fasta or gzipped fasta format
    - can be created with `harpy simulate {snpindel, inversion, cnv, translocation}`
    - [instructions here](simulate-variants.md)
- [optional] a file of 16-basepair barcodes to tag linked reads with
==- :icon-question: LRSIM differences
The original [LRSIM](https://github.com/aquaskyline/LRSIM) is a lengthy Perl script that, like Harpy, outsources
to various other programs (SURVIVOR, DWGSIM, samtools, msort) and acts as a workflow through these programs. The Harpy
version of LRSIM keeps only the novel LRSIM code that creates linked reads from reads simulated by DWGSIM. The
rest of LRSIM's components are reincorporated into the Snakemake workflow governing the `simulate linkedreads`
module, while removing the SURVIVOR part since `harpy simulate {snpindel,inversion,...}` are used for that purpose.
#### Notable differences
- dependencies are expected to be on the `PATH`, not hardcoded to the folder LRSIM is running from
- `-r` parameter changed to folder prefix since Harpy uses `-g` for the haplotypes
- outputs are coded a little differently for flexibility (and use the `-r` parameter for some parts)
- SURVIVOR variant simulation functionality removed entirely
- DWGSIM, samtools, msort, and extractReads functionality moved into Harpy workflow
- DWGSIM call does not add additional variants
- uses newer version of DWGSIM
===

You may want to benchmark haplotag data on different kinds of genomic variants. To
do that, you'll need *known* variants (like those created by `harpy simulate`) and
linked-read sequences. This module will create (diploid) linked-read sequences from two genome haplotypes.
To accomplish this, Harpy uses a modified version of [LRSIM](https://github.com/aquaskyline/LRSIM),
and converts the LRSIM 10X-style output into Haplotag-style reads. To simulate linked reads, use:

```bash usage
harpy simulate linkedreads OPTIONS... HAP1_GENOME HAP2_GENOME
```
```bash example
harpy simulate linkedreads -t 4 -n 2  -l 100 -p 50  data/genome.hap1.fasta data/genome.hap2.fasta
```

## :icon-terminal: Running Options
In addition to the [common runtime options](/commonoptions.md), the `simulate linkedreads` module is configured using these command-line arguments:

| argument       | short name | type        |    default    | required | description                                                                                     |
|:---------------|:----------:|:------------|:-------------:|:--------:|:------------------------------------------------------------------------------------------------|
| `HAP1_GENOME`       |            | file path |       | **yes**  | Haplotype 1 of the diploid genome to simulate reads   |
| `HAP2_GENOME`       |            | file path |       | **yes**  | Haplotype 1 of the diploid genome to simulate reads   |
| `--outer-distance`  |    `-d`    | integer   | 350   |   | Outer distance between paired-end reads (bp)                 |
| `--distance-sd`     |    `-i`    | integer   |  15   |        | Standard deviation of read-pair distance                |
| `--barcodes`        |    `-b`    | file path |  [10X barcodes](https://github.com/aquaskyline/LRSIM/blob/master/4M-with-alts-february-2016.txt)   |        | File of linked-read barcodes to add to reads   |
| `--read-pairs`      |    `-n`    | number    |  600  |   | Number of read pairs to simulate, in millions       |
| `--molecule-length` |    `-l`    | integer   |  100  |   | Mean molecule length (kbp)                          |
| `--patitions`       |    `-p`    | integer   |  1500 |   | How many partitions to generate (×1000)             |
| `--molecules-per`   |    `-m`    | integer   |   10  |   | Average number of molecules per partition           |

## Barcodes
Barcodes, if provided, must be given as 16-basepair nucleotide sequences, one per line. If not provided,
Harpy will download the standard 10X Genomics `4M-with-alts-february-2016.txt` barcode set from the LRSIM
repository and use those. The barcode file should look like:
``` input.barcodes.txt
ATATGTACTCATACCA
GGATGTACTCATTCCA
TCACGTACTCATACCA
etc...
```
Harpy will convert the simulated 10X-style reads, where the 16-basepair barcode is at the beginning of read 1,
to haplotag format, where the barcode is coded in the sequence header under the `BX:Z` tag with the format
`AxxCxxBxxDxx`, where `xx` is a number between `00` and `96`. Using this format, a `00` would invalidate the
entire barcode due to a segment failing to be associated with a beadtag segment. In the simulated data, since
10X barcodes don't feature segments, failure to associate the first 16 bases of read 1 with barcodes provided
to `--barcodes` will appear as `BX:Z:A00C00B00D00`. The original 10X barcode (or first 16 bases of read 1)
will be removed from the sequence and stored in the `TX:Z` sequence header tag, e.g. `TX:Z:ATATGTACTCATACCA`.
The paired reverse read will also have these tags.

## Choosing parameters
LRSIM does internal calculations to determine the number of reads per molecule based on `--read-pairs`,
`--partitions`, and `--molecules-per`. It would be helpful to see this equation so as to guide your
decisions for these parmaters. The equation is:

$$
ReadsPerMolecule = int( 0.499 + ( X \times 1,000,000 ) / ( T \times 1,000 / D ) / M / D )
$$