---
label: Variants
description: Simulate snps, indels, inversions, cnv, translocations
icon: dot
#visibility: hidden
order: 6
---

# :icon-flame: Simulate Genomic Variants
Simulate snps, indels, inversions, cnv, translocations

===  :icon-checklist: You will need
- a reference genome in fasta or gzipped fasta format
===

You may want to benchmark haplotag data on different kinds of genomic variants. To
do that, you'll need *known* variants, and typically simulations are how you achieve
that. This series of modules simulates genomic variants onto a genome, either randomly
or specific variants provided in VCF files. The simulator Harpy uses,
[simuG](https://github.com/yjx1217/simuG), can only simulate one type of
variant at a time and each variant type has their own set of parameters. This page
is divided by variant types to help you navigate the process. The general usage
for simulating variants is:

```bash usage
harpy simulate variant OPTIONS... INPUT_GENOME
```
```bash example
harpy simulate inversion -n 10 --min-size 1000 --max-size 50000  path/to/genome.fasta
```
There are 4 submodules with very obvious names:

| submodule | what it does |
|:----------|:-------------|
|`snpindel` | simulates single nucleotide polymorphisms (snps) and insertion-deletions (indels) |
| `inversion` | simulates inversions |
| `cnv` | simulates copy number variants |
| `translocation` | simulates translocations |

## :icon-terminal: Running Options
While there are serveral differences between the submodule command line options, each has available all the
[common runtime options](/commonoptions.md) like other Harpy modules. Each requires and input genome at the
end of the command line, and each requires either a `--count` of variants to randomly simulate, or a `--vcf` of
specific variants to simulate. There are also these unifying options among the different variant types:

| argument | short name | type | required | description |
| :-----|:-----|:-----|:-----|:-----|
| `INPUT_GENOME`           |            | file path  |   **yes**  | The haploid genome to simulate variants onto   |
| `--prefix` | | string | | Naming prefix for output files (default: `sim.{module_name}`)|
| `--exclude-chr` | `-e` | file path | | Text file of chromosomes to avoid, one per line |
| `--centromeres` | `-c` | file path | | GFF3 file of centromeres to avoid |
| `--genes` | `-g` | file path | | GFF3 file of genes to avoid simulating over (see `snpindel` for caveat) |
| `--heterozygosity` | `-z` | float between [0,1] | | [% heterozygosity to simulate diploid later](#simulate-diploids) (default: `0`) |
| `--randomseed` |  | integer > 0 |  | Random seed for simulation |

==- snpindel
### snpindel
The snp and indel variants are combined in this module because `simuG` allows simulating them together. The
ratio parameters control different things for snp and indel variants and have special meanings when setting
the value to either `9999` or `0` :
- `--titv-ratio`
    - `9999`: transitions only
    - `0`: transversions only
- `--indel-ratio`
    - `9999`: insertions only
    - `0`: deletions only

| argument          | short name | type       | default |  description                                                 |
|:------------------|:----------:|:-----------|:-------:|:-------------------------------------------------------------|
| `--snp-vcf`| `-v` | file path | | VCF file of known snps to simulate |
| `--indel-vcf` | `-i` | file path | | VCF file of known indels to simulate |
| `--snp-count` | `-n` | integer | 0 | Number of random snps to simluate |
| `--indel-count` |  `-m` | integer | 0 | Number of random indels to simluate |
| `--titv-ratio` | `-r` | float  | 0.5 | Transition/Transversion ratio for snps |
| `--indel-ratio` | `-d` | float  |  1 | Insertion/Deletion ratio for indels |
| `--indel-size-alpha` | `-a` | float |  2.0 | Exponent Alpha for power-law-fitted indel size distribution|
| `--indel-size-constant` | `-l` | float | 0.5 | Exponent constant for power-law-fitted indel size distribution |
| `--snp-gene-constraints` | `-y` | string | | How to constrain randomly simulated SNPs {`noncoding`,`coding`,`2d`,`4d`} when using `--genes`|

==- inversion
### inversion
Inversions are when a section of a chromosome appears in the reverse orientation.

| argument          | short name | type       | default |  description     |
|:------------------|:----------:|:-----------|:-------:|:----------------|
| `--vcf` | `-v` | file path |  |  VCF file of known inversions to simulate |
| `--count`| `-n` | integer | 0 |  Number of random inversions to simluate |
| `--min-size` | `-m` | integer | 1000 | Minimum inversion size (bp) |
| `--max-size` | `-x` | integer | 100000 | Maximum inversion size (bp) |

===
## Simulate known variants

## Simulate diploids
