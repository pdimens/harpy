---
label: Read Mapping
icon: quote
order: 5
---

# Mapping Reads onto a Reference Genome
You can map reads onto genome assemblies with Harpy by calling the `align` module:
```bash
harpy align OPTIONS...
```
|||  :icon-checklist: You will need
- at least 4 cores/threads available
- a genome assembly in FASTA format
- b/gzipped fastq sequence files
    - paired-end files
    - file extension is either `.fastq.gz` or `.fq.gz` (do not mix)
    - forward-reverse is noted as either `.R1.`/`.R2.` or `_R1.`/`_R2.` (do no mix)
        - _e.g._ `samplename.R1.fq.gz` and `samplename.R2.fq.gz`
        - _e.g._ `samplename_R1.fq.gz` and `samplename_R2.fq.gz`
        - or the same but ending with `.fastq.gz`, but don't mix and match
|||

## Running Options
| argument         | short name | type        | default | required | description                                    |
|:-----------------|:----------:|:------------|:-------:|:--------:|:-----------------------------------------------|
| `--genome`       |    `-g`    | file path   |         | **yes**  | Genome assembly for read mapping               |
| `--dir`          |    `-d`    | folder path |         | **yes**  | Directory with sample sequences                |
| `--ema-bins`     |    `-e`    | integer     |   500   |    no    | Number of barcode bins for EMA                 |
| `--bwa`          |    `-b`    | toggle      |         |    no    | Use BWA MEM instead of EMA                     |
| `--extra-params` |    `-x`    | string      |         |    no    | Additional EMA-align/BWA parameters, in quotes |
| `--threads`      |    `-t`    | integer     |    4    |    no    | Number of threads to use                       |
| `--snakemake`    |    `-s`    | string      |         |    no    | Additional Snakemake options, in quotes        |
| `--help`         |            |             |         |          | Show the module docstring                      |

----

## Workflows
### EMA
||| EMA details
- **recommended**, [see below](readmapping.md/#why-ema)
- leverages the BX barcode information to improve mapping
- slower
- lots of temporary files
|||

Since [EMA](https://github.com/arshajii/ema) does extra things to account for barcode information, the EMA workflow is a bit more complicated under the hood. Reads with barcodes are aligned using EMA and reads without valid barcodes are separately mapped using BWA before merging all the alignments together again. EMA will mark duplicates within alignments, but the BWA alignments need duplicates marked manually using [sambamba](https://lomereiter.github.io/sambamba/). Thankfully, you shouldn't need to worry about any of these details.

```mermaid
graph LR
    A((count beadtags)) --> B((EMA preprocess))
    B-->C((EMA align barcoded))
    C-->D((sort alignments))
    D-->E((merge alignments))
    E-->F((merge alignemnts))
    IDX((index genome))-->C
    IDX-->Z((BWA align unbarcoded))
    Z-->Y((sort alignments))
    Y-->X((mark duplicates))
    X-->F
    F-->K((sort alignments))
    K-->J((alignment reports))
    K-->G((convert to BED))
    G-->H((calculate reads per BX))
    G-->L((calculate BX size))
    G-->I((calculate genomic coverage))
```

### BWA
||| BWA details
- ignores barcode information
- might be preferred depending on experimental design
- faster
- no temporary files
|||

The [BWA MEM](https://github.com/lh3/bwa) workflow is substantially simpler than the EMA workflow and maps all reads against the reference genome, no muss no fuss. Duplicates are marked at the end using [sambamba](https://lomereiter.github.io/sambamba/).

```mermaid
graph LR
    A((index genome)) --> B((align to genome))
    B-->C((sort alignments))
    C-->D((mark duplicates))
    D-->E((alignment reports))
    D-->F((convert to BED))
    F-->G((calculate genomic coverage))
```

## Why EMA?
The original haplotag manuscript uses BWA to map reads, but the authors have since then recommended the use of EMA (EMerald Aligner) for most applications. EMA is barcode-aware, meaning it considers sequences with identical barcodes to have originated from the same molecule, and therefore has higher mapping accuracy than using BWA. Here's a comparison from the [EMA manuscript](https://www.biorxiv.org/content/10.1101/220236v1):
![EMA Publication figure 3](/static/EMA.fig3.png)