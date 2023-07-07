---
label: Demultiplex
description: Demultiplex raw sequences into haplotag barcoded samples
icon: versions
order: 7
---

# :icon-versions: Demultiplex Raw Sequences

===  :icon-checklist: You will need
- at least 2 cores/threads available
- gzipped forward/reverse reads from an Illumina sequencer
===


You may want to phase your genotypes into haplotypes, as haplotypes tend to be more informative
than unphased genotypes (higher polymorphism, captures relationship between genotypes). Phasing
genotypes into haplotypes requires alignment files, such as those produced by `harpy align` and
a variant call file, such as those produced by `harpy variants` or `harpy impute`. **Phasing only
works on SNP data**, and will not work for structural variants produced by `LEVIATHAN`. You can phase genotypes into haplotypes with Harpy using the `phase` module:

```bash usage
harpy demultiplex OPTIONS... 
```
```bash example
harpy phase --threads 20 --directory Variants/variants.raw.bcf
```


## :icon-terminal: Running Options
| argument              | short name | type            | default | required | description                                                                                     |
|:----------------------|:----------:|:----------------|:-------:|:--------:|:------------------------------------------------------------------------------------------------|
| `--vcf`               |    `-v`    | file path       |         | **yes**  | Path to BCF/VCF file                                                                            |
| `--directory`               |    `-d`    | folder path     |         | **yes**  | Directory with sequence alignments                                                              |
| `--molecule-distance` |    `-m`    | integer         |  20000  |    no    | Base-pair distance dilineating separate molecules                                               |
| `--prune-threshold`   |    `-p`    | integer (0-100) |    7    |    no    | PHRED-scale (%) threshold for pruning low-confidence SNPs                                       |
| `--extra-params`      |    `-x`    | string          |         |    no    | Additional Hapcut2 arguments, in quotes                                                         |
| `--threads`           |    `-t`    | integer         |    4    |    no    | Number of threads to use                                                                        |
| `--snakemake`         |    `-s`    | string          |         |    no    | Additional [Snakemake](../snakemake/#adding-snakamake-parameters) options, in quotes |
| `--quiet`             |    `-q`    | toggle          |         |    no    | Supressing Snakemake printing to console                                                        |
| `--help`              |            |                 |         |          | Show the module docstring                                                                       |


The molecule distance is and pruning thresholds are considered the most impactful parameters
for running HapCut2, therefore they are directly configurable from the command. The molecule distance
refers to the base-pair distance dilineating separate molecules. Feel free to play around with this number 
if you do not know the distance, as it's not clear how impactful this can be on the results. The pruning 
threshold refers to a PHRED-scale value between 0-1 (a percentage) for removing low-confidence SNPs from consideration. 
With Harpy, you configure this value as an integer between 0-100, which gets converted to a floating point
value between 0-1 internally (_i.e._ `-p 7` is equivalent to 0.07 threshold).

---
## :icon-git-pull-request: Phasing Workflow
+++ :icon-git-merge: details
Phasing is performed using [HapCut2](https://github.com/vibansal/HapCUT2). Most of the tasks cannot
be parallelized, but HapCut2 operates on a per-sample basis, so the workflow is parallelized
across all of your samples to speed things along.

```mermaid
graph LR
    A([split samples]) --> B([extractHAIRS])
    B-->C([LinkFragments])
    Z([sample alignments]) --> B
    Z-->C
    C-->D([phase blocks])
    B-->D
    A-->D
    D-->E([annotate BCFs])
    E-->F([index annotations])
    F-->G([merge annotations])
    E-->G
    A-->G
    D-->G
    G-->H([index merged annotations])
    H-->I([merge phased samples])
```

+++ :icon-file-directory: phasing output
The `harpy phase` module creates an `Phasing` directory with the folder structure below. `Sample1` is a generic sample name for demonstration purposes. Harpy will also write a record of the relevant
runtime parameters in `logs/phase.params`.

```
Phase/
├── variants.phased.bcf
├── variants.phased.bcf.csi
├── annotations
│   ├── Sample1.annot.gz
│   └── Sample1.annot.gz.tbi
├── annotations_merge
│   ├── Sample1.phased.annot.bcf
│   └── Sample1.phased.annot.bcf.csi
├── extractHairs
│   ├── Sample1.unlinked.frags
│   └── logs
│       └── Sample1.unlinked.log
├── input
│   ├── header.names
│   ├── Sample1.bcf
│   └── Sample1.het.bcf
├── linkFragments
│   ├── Sample1.linked.frags
│   └── logs
│       └── Sample1.linked.log
└── phaseBlocks
    ├── Sample1.blocks
    ├── Sample1.blocks.phased.VCF
    └── logs
        └── Sample1.blocks.phased.log

```

| item | description |
|:---|:---|
| `variants.phased.bcf*` | final vcf output of HapCut2 with all samples merged into a single file (with .csi index) |
| `annotations/` | phased vcf annotated with phased blocks |
| `annotations_merge/` | merged vcf of annotated and original vcf |
| `extractHairs/` | output from `extractHairs` |
| `extractHairs/logs/` | everything HapCut2's `extractHairs` prints to `stderr` |
| `input/head.names` | extra file harpy creates to support new INFO fields in the phased VCF |
| `input/*.bcf` | vcf of a single sample from the original multi-sample input vcf |
| `input/*.het.bcf` | vcf of heterozygous loci of a single sample from the original multi-sample input vcf |
| `linkFragments/` | results from HapCut2's `linkFragments` |
| `linkFragments/logs` | everything `linkFragments` prints to `stderr` |
| `phaseBlocks/*.blocks*` | output from HapCut2 |
| `phaseBlocks/logs` | everything HapCut2 prints to `stderr` |

+++ :icon-code-square: HapCut2 parameters
By default, Harpy runs `HAPCUT2` with these parameters (excluding inputs and outputs):
```bash
HAPCUT2 --nf 1 --error_analysis_mode 1 --call_homozygous 1 --outvcf 1
```
Below is a list of all `HAPCUT2` command line options, excluding those Harpy already uses or those made redundant by Harpy's implementation of HapCut2.
These are taken directly from running `HAPCUT2 --help`.

``` hapcut2 arguments
Haplotype Post-Processing Options:
--skip_prune, --sp <0/1>:           skip default likelihood pruning step (prune SNPs after the fact using column 11 of the output). default: 0
--discrete_pruning, --dp <0/1>:     use discrete heuristic to prune SNPs. default: 0
--error_analysis_mode, --ea <0/1>:  compute switch confidence scores and print to haplotype file but don't split blocks or prune. default: 0

Advanced Options:
--max_iter, --mi <int> :            maximum number of global iterations. Preferable to tweak --converge option instead. default: 10000
--maxcut_iter, --mc <int> :         maximum number of max-likelihood-cut iterations. Preferable to tweak --converge option instead. default: 10000
```
+++
