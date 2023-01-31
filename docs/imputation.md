# Impute Genotypes from Genotype Likelihoods
You can impute genotypes with Harpy by calling the `impute` module:
```bash
harpy impute OPTIONS...
```
You will need:
- a tab-delimited parameter file 
    - create with `harpy impute --init`
    - modify the file with parameters suitable for your study
- a variant call format file with genotype likelihoods
    - accepted formats: `.vcf`, `.vcf.gz`, `.bcf`
- sequence alignments, in `.bam` format

## Running Options
| argument       | short name | type        |    default    | required | description                                                   |
|:---------------|:----------:|:------------|:-------------:|:--------:|:--------------------------------------------------------------|
| `--init`       |    `-i`    | toggle      |               |          | Create example parameter file and exit                        |
| `--vcf`        |    `-v`    | file path   |               |   **yes**    | Path to VCF/BCF file                                          |
| `--directory`  |    `-d`    | folder path |               |   **yes**    | Directory with sequence alignments                            |
| `--parameters` |    `-p`    | file path   | stitch.params |   **yes**    | STITCH parameter file (tab-delimited)                         |
| `--filter`     |    `-f`    | toggle      |               |    no    | Filter `--vcf` file to keep SNPs with Quality>20 and Depth>10 |
| `--threads`    |    `-t`    | integer     |       4       |    no    | Number of threads to use                                      |
| `--snakemake`  |    `-s`    | string      |               |    no    | Additional Snakemake options, in quotes                       |
| `--help`       |            |             |               |          | Show the module docstring                                     |

## Parameter file
Typically, one runs STITCH multiple times, exploring how results vary with
different model parameters. The solution Harpy uses for this is to have the user
provide a tab-delimited dataframe file where the columns are the 5 STITCH model 
parameters and the rows are the values for those parameters. The parameter file 
is required and can be created manually or with `harpy impute --init`.
If created using harpy, the resulting file includes largely meaningless values 
that you will need to adjust for your study. The parameter must follow a particular format:
- tab or comma delimited
- column order doesn't matter, but all 5 column names must be present
- header row present with the specific column names below
    - all column names begin with a lowercase character
| column name |  value type  |             accepted values             | description                                                           |
|:------------|:------------:|:---------------------------------------:|:----------------------------------------------------------------------|
| model       |     text     | pseudoHaploid, diploid, diploid-inbred  | The STITCH model/method to use                                        |
| useBX       | text/boolean | true, false, yes, no (case insensitive) | Whether to incorporate beadtag information                            |
| k           |   integer    |                   ≥ 1                   | Number of founder haplotypes                                          |
| s           |   integer    |                   ≥ 1                   | Number of instances of the founder haplotypes to average results over |
| nGen        |   integer    |                   ≥ 1                   | Estimated number of generations since founding                        |

### example file
<!-- tabs:start -->

#### **tab-delimited**

This file is tab-delimited, note the column names:

```
model   useBX   k       s       nGen
pseudoHaploid   TRUE    10      5       50
pseudoHaploid   TRUE    10      1       50
pseudoHaploid   TRUE    15      10      100
```

#### **table-view**

This is the table view of the tab-delimited file, shown here for clarity.

| model         | useBX | k  | s  | nGen |
|:--------------|:------|:---|:---|:-----|
| pseudoHaploid | TRUE  | 10 | 5  | 50   |
| pseudoHaploid | TRUE  | 10 | 1  | 50   |
| pseudoHaploid | TRUE  | 15 | 10 | 100  |

<!-- tabs:end -->

## STITCH Workflow
[STITCH](https://github.com/rwdavies/STITCH) is a genotype imputation software developed for use in
the R programming language. It has quite a few model parameters that can be tweaked, but HARPY only
focuses on a small handful that have the largest impact on the quality of the results. Imputation is
performed on a per-contig (or chromosome) level, so Harpy automatically iterates over the contigs
present in the input variant call file. Using the magic of Snakemake, Harpy will automatically
iterate over these model parameters.

```mermaid
graph LR
    A((count contigs)) --> B((split contigs))
    B-->C((keep biallelic SNPs))
    C-->D((convert to STITCH format))
    D-->E((STITCH imputation))
    E-->F((merge output))
    G((create file list))-->E
```