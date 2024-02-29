---
label: BWA
description: Align haplotagged sequences with BWA MEM
icon: dot
order: 5
---

# :icon-quote: Map Reads onto a genome with BWA MEM
===  :icon-checklist: You will need
- at least 4 cores/threads available
- a genome assembly in FASTA format
- paired-end fastq sequence file with the [proper naming convention](/haplotagdata/#naming-conventions) (gzipped recommended)
===

Once sequences have been trimmed and passed through other QC filters, they will need to
be aligned to a reference genome. This module within Harpy expects filtered reads as input,
such as those derived using `harpy qc`. You can map reads onto a genome assembly with Harpy 
using the `align` module:

```bash usage
harpy align bwa OPTIONS... INPUTS...
```
```bash example
harpy align bwa --genome genome.fasta Sequences/ 
```

## :icon-terminal: Running Options
In addition to the [common runtime options](/commonoptions.md), the `harpy align bwa` module is configured using these command-line arguments:

| argument           | short name | type                  | default | required | description                                           |
|:-------------------|:----------:|:----------------------|:-------:|:--------:|:------------------------------------------------------|
| `INPUTS`           |            | file/directory paths  |         | **yes**  | Files or directories containing [input FASTQ files](/commonoptions.md#input-arguments)     |
| `--genome`         |    `-g`    | file path             |         | **yes**  | Genome assembly for read mapping                      |
| `--molecule-distance` |    `-m`    | integer         |  100000  |    no    | Base-pair distance threshold to separate molecules      |
| `--quality-filter` |    `-f`    | integer (0-40)        |   30    |    no    | Minimum `MQ` (SAM mapping quality) to pass filtering  |
| `--method`         |    `-m`    | choice [`bwa`, `ema`] |   bwa   |    no    | Which aligning software to use                        |
| `--extra-params`   |    `-x`    | string                |         |    no    | Additional EMA-align/BWA arguments, in quotes         |

### Molecule distance
The `--molecule-distance` option is used during the BWA alignment workflow
to assign alignments a unique Molecular Identifier `MI:i` tag based on their
 haplotag barcode and the distance threshold you specify. See 
[haplotag data](/haplotagdata/#barcode-thresholds) for more information on
what this value does. 

## Quality filtering
The `--quality` argument filters out alignments below a given $MQ$ threshold. The default, `30`, keeps alignments
that are at least 99.9% likely correctly mapped. Set this value to `1` if you only want alignments removed with
$MQ = 0$ (0% likely correct). You may also set it to `0` to keep all alignments for diagnostic purposes.
The plot below shows the relationship between $MQ$ score and the likelihood the alignment is correct and will serve to help you decide
on a value you may want to use. It is common to remove alignments with $MQ <30$ (<99.9% chance correct) or $MQ <40$ (<99.99% chance correct).

==- What is the $MQ$ score?
Every alignment in a BAM file has an associated mapping quality score ($MQ$) that informs you of the likelihood 
that the alignment is accurate. This score can range from 0-40, where higher numbers mean the alignment is more
likely correct. The math governing the $MQ$ score actually calculates the percent chance the alignment is ***incorrect***: 
$$
\%\ chance\ incorrect = 10^\frac{-MQ}{10} \times 100\\
\text{where }0\le MQ\le 40
$$
You can simply subtract it from 100 to determine the percent chance the alignment is ***correct***:
$$
\%\ chance\ correct = 100 - \%\ chance\ incorrect\\
\text{or} \\
\%\ chance\ correct = (1 - 10^\frac{-MQ}{10}) \times 100
$$

[!embed el="embed"](//plotly.com/~pdimens/7.embed)
===

## Marking PCR duplicates
Harpy uses `samtools markdup` to mark putative PCR duplicates. By using the `--barcode-tag BX`
option, it considers the linked-read barcode for more accurate duplicate detection. Duplicate
marking also uses the `-S` option to mark supplementary (chimeric) alignments as duplicates
if the primary alignment was marked as a duplicate. Duplicates get marked but **are not removed**.

----

## :icon-git-pull-request: BWA workflow
+++ :icon-git-merge: details
- default aligner
- ignores (but retains) barcode information
- fast

The [BWA MEM](https://github.com/lh3/bwa) workflow is much simpler and faster than the EMA workflow
and maps all reads against the reference genome. Duplicates are marked using `samtools markdup`.
The `BX:Z` tags in the read headers are still added to the alignment headers, even though barcodes
are not used to inform mapping. The `-m` threshold is used for alignment molecule assignment.

```mermaid
graph LR
    Z([trimmed reads]) --> B
    A([index genome]) --> B([align to genome])
    B-->C([sort alignments])
    C-->D([mark duplicates])
    D-->E([assign molecules])
    E-->F([alignment metrics])
    D-->G([barcode stats])
    G-->F
    subgraph markdp [mark duplicates via `samtools`]
        direction LR
        collate-->fixmate
        fixmate-->sort
        sort-->markdup
    end
```
+++ :icon-file-directory: BWA output
The `harpy align` module creates an `Align/bwa` directory with the folder structure below. `Sample1` is a generic sample name for demonstration purposes.
The resulting folder also includes a `workflow` directory (not shown) with workflow-relevant runtime files and information.
```
Align/bwa
├── Sample1.bam
├── Sample1.bam.bai
├── logs
│   └── markduplicates
│       └── Sample1.markdup.log
└── reports
    ├── bwa.stats.html
    ├── BXstats
    │   ├── Sample1.bxstats.html
    │   └── data
    │       └── Sample1.bxstats.gz
    └── coverage
        ├── Sample1.gencov.html
        └── data
            └── Sample1.gencov.gz


```

| item     | description                                                                                                 |
|:---------|:------------------------------------------------------------------------------------------------------------|
| `*.bam`                             | sequence alignments for each sample                                              |
| `*.bai`                             | sequence alignment indexes for each sample                                       |
| `logs/markduplicates`               | stats provided by `samtools markdup`                                             |
| `reports/`                          | various counts/statistics/reports relating to sequence alignment                 |
| `reports/bwa.stats.html`            | report summarizing `samtools flagstat and stats` results across all samples from `multiqc` |
| `reports/reads.bxstats.html`        | interactive html report summarizing valid vs invalid barcodes across all samples | 
| `reports/BXstats/*.bxstats.html`    | interactive html report summarizing inferred molecule size                       | 
| `reports/coverage/*.html`           | summary plots of alignment coverage per contig                                   |
| `reports/coverage/data/*.gencov.gz` | output from samtools bedcov from all alignments, used for plots                  |
| `reports/BXstats/`                  | reports summarizing molecule size and reads per molecule                         |
| `reports/BXstats/data/`             | tabular data containing the information used to generate the BXstats reports     |

+++ :icon-code-square: BWA parameters
By default, Harpy runs `bwa` with these parameters (excluding inputs and outputs):
```bash
bwa mem -C -R "@RG\tID:samplename\tSM:samplename"
```

Below is a list of all `bwa mem` command line arguments, excluding those Harpy already uses or those made redundant by Harpy's implementation of BWA.
These are taken directly from the [BWA documentation](https://bio-bwa.sourceforge.net/bwa.shtml).
```bwa arguments
-k INT 	Minimum seed length. Matches shorter than INT will be missed. The alignment speed is usually insensitive to this value unless it significantly deviates 20. [19]
-w INT 	Band width. Essentially, gaps longer than INT will not be found. Note that the maximum gap length is also affected by the scoring matrix and the hit length, not solely determined by this option. [100]
-d INT 	Off-diagonal X-dropoff (Z-dropoff). Stop extension when the difference between the best and the current extension score is above |i-j|*A+INT, where i and j are the current positions of the query and reference, respectively, and A is the matching score. Z-dropoff is similar to BLAST’s X-dropoff except that it doesn’t penalize gaps in one of the sequences in the alignment. Z-dropoff not only avoids unnecessary extension, but also reduces poor alignments inside a long good alignment. [100]
-r FLOAT 	Trigger re-seeding for a MEM longer than minSeedLen*FLOAT. This is a key heuristic parameter for tuning the performance. Larger value yields fewer seeds, which leads to faster alignment speed but lower accuracy. [1.5]
-c INT 	Discard a MEM if it has more than INT occurence in the genome. This is an insensitive parameter. [10000]
-P 	In the paired-end mode, perform SW to rescue missing hits only but do not try to find hits that fit a proper pair.
-A INT 	Matching score. [1]
-B INT 	Mismatch penalty. The sequence error rate is approximately: {.75 * exp[-log(4) * B/A]}. [4]
-O INT 	Gap open penalty. [6]
-E INT 	Gap extension penalty. A gap of length k costs O + k*E (i.e. -O is for opening a zero-length gap). [1]
-L INT 	Clipping penalty. When performing SW extension, BWA-MEM keeps track of the best score reaching the end of query. If this score is larger than the best SW score minus the clipping penalty, clipping will not be applied. Note that in this case, the SAM AS tag reports the best SW score; clipping penalty is not deducted. [5]
-U INT 	Penalty for an unpaired read pair. BWA-MEM scores an unpaired read pair as scoreRead1+scoreRead2-INT and scores a paired as scoreRead1+scoreRead2-insertPenalty. It compares these two scores to determine whether we should force pairing. [9]
-T INT 	Don’t output alignment with score lower than INT. This option only affects output. [30]
-a 	Output all found alignments for single-end or unpaired paired-end reads. These alignments will be flagged as secondary alignments.
-H 	Use hard clipping ’H’ in the SAM output. This option may dramatically reduce the redundancy of output when mapping long contig or BAC sequences.
```
+++ :icon-graph: reports
These are the summary reports Harpy generates for this workflow. You may right-click
the images and open them in a new tab if you wish to see the examples in better detail.
||| Depth and coverage
Reports the depth of alignments in 10kb windows.
![reports/coverage/*.html](/static/report_align_coverage.png)
||| BX validation
Reports the number of valid/invalid barcodes in the alignments.
![reports/reads.bxstats.html](/static/report_align_bxstats.png)
||| Molecule size
Reports the inferred molecule sized based on barcodes in the alignments.
![reports/BXstats/*.bxstats.html](/static/report_align_bxmol.png)
||| Alignment stats
Reports the general statistics computed by samtools `stats` and `flagstat`
![reports/samtools_*stat/*html](/static/report_align_flagstat.png)
|||

+++
