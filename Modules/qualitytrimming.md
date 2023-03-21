---
label: Quality Trimming
description: Quality trim haplotagged sequences with Harpy
icon: codescan-checkmark
order: 6
---

# Quality Trimming Sequence Data
===  :icon-checklist: You will need
- at least 2 cores/threads available
- paired-end b/gzipped fastq sequence files
===

Raw sequences are not suitable for downstream analyses. They have sequencing adapters,
index sequences, regions of poor quality, etc. The first step of any genetic sequence
analyses is to remove these adapters and trim poor quality data. You can remove adapters
and quality trim sequences using the `trim` module:

```bash usage
harpy trim OPTIONS... 
```

```bash example
harpy trim --dir Sequences_Raw/ --threads 20 
```

## Running Options
| argument         | short name | type        | default | required | description                                                            |
|:-----------------|:----------:|:------------|:-------:|:--------:|:-----------------------------------------------------------------------|
| `--dir`          |    `-d`    | folder path |         | **yes**  | Directory with sequence alignments                                     |
| `--max-length`   |    `-l`    | integer     |   150   |    no    | Maximum length to trim sequences down to                               |
| `--extra-params` |    `-x`    | string      |         |    no    | Additional fastp arguments, in quotes                                 |
| `--threads`      |    `-t`    | integer     |    4    |    no    | Number of threads to use                                               |
| `--snakemake`    |    `-s`    | string      |         |    no    | Additional Snakemake options, in quotes ([more info](../getstarted.md/#adding-additional-snakamake-parameters)) |
| `--help`         |            |             |         |          | Show the module docstring                                              |

==- fastp arguments
Below is a list of all `fastp` command line options, excluding those Harpy already uses or those made redundant by Harpy's implementation of fastp.
These are taken directly from the [fastp documentation](https://github.com/OpenGene/fastp).
``` fastp arguments
  # I/O options
      --failed_out                     specify the file to store reads that cannot pass the filters. (string [=])
      --overlapped_out                 for each read pair, output the overlapped region if it has no any mismatched base. (string [=])
  -m, --merge                          for paired-end input, merge each pair of reads into a single read if they are overlapped. The merged reads will be written to the file given by --merged_out, the unmerged reads will be written to the files specified by --out1 and --out2. The merging mode is disabled by default.
      --merged_out                     in the merging mode, specify the file name to store merged output, or specify --stdout to stream the merged output (string [=])
      --include_unmerged               in the merging mode, write the unmerged or unpaired reads to the file specified by --merge. Disabled by default.
  -6, --phred64                        indicate the input is using phred64 scoring (it'll be converted to phred33, so the output will still be phred33)
  -z, --compression                    compression level for gzip output (1 ~ 9). 1 is fastest, 9 is smallest, default is 4. (int [=4])
      --reads_to_process               specify how many reads/pairs to be processed. Default 0 means process all reads. (int [=0])
      --fix_mgi_id                     the MGI FASTQ ID format is not compatible with many BAM operation tools, enable this option to fix it.
  
  # adapter trimming options
  -A, --disable_adapter_trimming       adapter trimming is enabled by default. If this option is specified, adapter trimming is disabled
  -a, --adapter_sequence               the adapter for read1. For SE data, if not specified, the adapter will be auto-detected. For PE data, this is used if R1/R2 are found not overlapped. (string [=auto])
      --adapter_sequence_r2            the adapter for read2 (PE data only). This is used if R1/R2 are found not overlapped. If not specified, it will be the same as <adapter_sequence> (string [=])
      --adapter_fasta                  specify a FASTA file to trim both read1 and read2 (if PE) by all the sequences in this FASTA file (string [=])
    
  # global trimming options
  -f, --trim_front1                    trimming how many bases in front for read1, default is 0 (int [=0])
  -t, --trim_tail1                     trimming how many bases in tail for read1, default is 0 (int [=0])
  -F, --trim_front2                    trimming how many bases in front for read2. If it's not specified, it will follow read1's settings (int [=0])
  -T, --trim_tail2                     trimming how many bases in tail for read2. If it's not specified, it will follow read1's settings (int [=0])

  # duplication evaluation and deduplication
  -D, --dedup                          enable deduplication to drop the duplicated reads/pairs
      --dup_calc_accuracy              accuracy level to calculate duplication (1~6), higher level uses more memory (1G, 2G, 4G, 8G, 16G, 24G). Default 1 for no-dedup mode, and 3 for dedup mode. (int [=0])
      --dont_eval_duplication          don't evaluate duplication rate to save time and use less memory.

  # polyG tail trimming, useful for NextSeq/NovaSeq data
      --poly_g_min_len                 the minimum length to detect polyG in the read tail. 10 by default. (int [=10])
  -G, --disable_trim_poly_g            disable polyG tail trimming, by default trimming is automatically enabled for Illumina NextSeq/NovaSeq data

  # polyX tail trimming
  -x, --trim_poly_x                    enable polyX trimming in 3' ends.
      --poly_x_min_len                 the minimum length to detect polyX in the read tail. 10 by default. (int [=10])
  
  # per read cutting by quality options
  -5, --cut_front                      move a sliding window from front (5') to tail, drop the bases in the window if its mean quality < threshold, stop otherwise.
  -3, --cut_tail                       move a sliding window from tail (3') to front, drop the bases in the window if its mean quality < threshold, stop otherwise.
  -W, --cut_window_size                the window size option shared by cut_front, cut_tail or cut_sliding. Range: 1~1000, default: 4 (int [=4])
  -M, --cut_mean_quality               the mean quality requirement option shared by cut_front, cut_tail or cut_sliding. Range: 1~36 default: 20 (Q20) (int [=20])
      --cut_front_window_size          the window size option of cut_front, default to cut_window_size if not specified (int [=4])
      --cut_front_mean_quality         the mean quality requirement option for cut_front, default to cut_mean_quality if not specified (int [=20])
      --cut_tail_window_size           the window size option of cut_tail, default to cut_window_size if not specified (int [=4])
      --cut_tail_mean_quality          the mean quality requirement option for cut_tail, default to cut_mean_quality if not specified (int [=20])
      --cut_right_window_size          the window size option of cut_right, default to cut_window_size if not specified (int [=4])
      --cut_right_mean_quality         the mean quality requirement option for cut_right, default to cut_mean_quality if not specified (int [=20])
  
  # quality filtering options
  -Q, --disable_quality_filtering      quality filtering is enabled by default. If this option is specified, quality filtering is disabled
  -q, --qualified_quality_phred        the quality value that a base is qualified. Default 15 means phred quality >=Q15 is qualified. (int [=15])
  -u, --unqualified_percent_limit      how many percents of bases are allowed to be unqualified (0~100). Default 40 means 40% (int [=40])
  -n, --n_base_limit                   if one read's number of N base is >n_base_limit, then this read/pair is discarded. Default is 5 (int [=5])
  -e, --average_qual                   if one read's average quality score <avg_qual, then this read/pair is discarded. Default 0 means no requirement (int [=0])

  
  # length filtering options
  -L, --disable_length_filtering       length filtering is enabled by default. If this option is specified, length filtering is disabled
  -l, --length_required                reads shorter than length_required will be discarded, default is 15. (int [=15])
      --length_limit                   reads longer than length_limit will be discarded, default 0 means no limitation. (int [=0])

  # low complexity filtering
  -y, --low_complexity_filter          enable low complexity filter. The complexity is defined as the percentage of base that is different from its next base (base[i] != base[i+1]).
  -Y, --complexity_threshold           the threshold for low complexity filter (0~100). Default is 30, which means 30% complexity is required. (int [=30])

  # filter reads with unwanted indexes (to remove possible contamination)
      --filter_by_index1               specify a file contains a list of barcodes of index1 to be filtered out, one barcode per line (string [=])
      --filter_by_index2               specify a file contains a list of barcodes of index2 to be filtered out, one barcode per line (string [=])
      --filter_by_index_threshold      the allowed difference of index barcode for index filtering, default 0 means completely identical. (int [=0])

  # base correction by overlap analysis options
  -c, --correction                     enable base correction in overlapped regions (only for PE data), default is disabled
      --overlap_len_require            the minimum length to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 30 by default. (int [=30])
      --overlap_diff_limit             the maximum number of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. 5 by default. (int [=5])
      --overlap_diff_percent_limit     the maximum percentage of mismatched bases to detect overlapped region of PE reads. This will affect overlap analysis based PE merge, adapter trimming and correction. Default 20 means 20%. (int [=20])

  # UMI processing
  -U, --umi                            enable unique molecular identifier (UMI) preprocessing
      --umi_loc                        specify the location of UMI, can be (index1/index2/read1/read2/per_index/per_read, default is none (string [=])
      --umi_len                        if the UMI is in read1/read2, its length should be provided (int [=0])
      --umi_prefix                     if specified, an underline will be used to connect prefix and UMI (i.e. prefix=UMI, UMI=AATTCG, final=UMI_AATTCG). No prefix by default (string [=])
      --umi_skip                       if the UMI is in read1/read2, fastp can skip several bases following UMI, default is 0 (int [=0])

  # overrepresented sequence analysis
  -p, --overrepresentation_analysis    enable overrepresented sequence analysis.
  -P, --overrepresentation_sampling    One in (--overrepresentation_sampling) reads will be computed for overrepresentation analysis (1~10000), smaller is slower, default is 20. (int [=20])

  # output splitting options
  -s, --split                          split output by limiting total split file number with this option (2~999), a sequential number prefix will be added to output name ( 0001.out.fq, 0002.out.fq...), disabled by default (int [=0])
  -S, --split_by_lines                 split output by limiting lines of each file with this option(>=1000), a sequential number prefix will be added to output name ( 0001.out.fq, 0002.out.fq...), disabled by default (long [=0])
  -d, --split_prefix_digits            the digits for the sequential number padding (1~10), default is 4, so the filename will be padded as 0001.xxx, 0 to disable padding (int [=4])
  
  # help
  -?, --help                           print this message
```
===

## FASTQ file format
There are a handful of "accepted" naming schemes for fastq file extensions, but Harpy only accepts a limited number of them, shown below.
The fastq files **must** be bgzipped or gzipped and be **consistent** with regards to the extensions and read-pair naming styles.
That is, you must only use `.fastq.gz` or only use `.fq.gz` for all files, and the same for `.1.`/`.2.` or `.F.`/`.R.` (adhere to a single row in the table below).
Notice that the read pair part differs from the [accepted fastq formats](readmapping.md/#fastq-file-format) for aligning reads.
=== acceptable formats

| forward-reverse notation | extension  | example forward          | example reverse         |
|:-------------------------|:-----------|:-------------------------|:------------------------|
| `.F` / `.R`                | `fastq.gz` | ` samplename.F.fastq.gz` | `samplename.R.fastq.gz` |
| `.F` / `.R`                | `fq.gz`    | `samplename.F.fq.gz`     | `samplename.R.fq.gz`    |
| `.1` / `.2`                | `fastq.gz` | `samplename.2.fastq.gz`  | `samplename.2.fastq.gz` |
| `.1` / `.2`                | `fq.gz`    | `samplename.1.fq.gz`     | `samplename.2.fq.gz`    |

===

---
## Trimming Workflow
+++ description
[Fastp](https://github.com/OpenGene/fastp) is an ultra-fast all-in-one adapter remover, deduplicator, 
and quality trimmer. Harpy uses it to remove adapters, low-quality bases, and trim sequences down to a particular
length (default 150bp). Harpy uses the fastp overlap analysis to identify adapters for removal and a sliding window
approach (`--cut-right`) to identify low quality bases. The workflow is quite simple.

```mermaid
graph LR
    A([fastp trim]) --> B([create reports])
```

+++ trimming output
The `harpy trim` module creates a `Trimming` directory with the folder structure below. `Sample1` and `Sample2` are generic sample names for demonstration purposes.
```
Trimming/
├── Sample1.R1.fq.gz
├── Sample1.R2.fq.gz
├── Sample2.R1.fq.gz
├── Sample2.R2.fq.gz
└── logs
    ├── err
    │   ├── Sample1.log
    │   └── Sample2.log
    ├── html
    │   ├── Sample1.html
    │   └── Sample2.html
    ├── json
    │   ├── Sample1.fastp.json
    │   └── Sample2.fastp.json
    └── trim.report.html
```
| item                    | description                                                              |
|:------------------------|:-------------------------------------------------------------------------|
| `*.R1.fq.gz`            | quality trimmed forward reads of the samples                             |
| `*.R1.fq.gz`            | quality trimmed reverse reads of the samples                             |
| `logs/`                 | all debug/diagnostic files that aren't the trimmed reads `fastp` creates |
| `logs/err`              | what fastp prints to `stderr` when running                               |
| `logs/html`             | interactive html reports `fastp` creates from quality trimming           |
| `logs/json`             | json representation of the data `fastp` uses to create the html reports  |
| `logs/trim.report.html` | a report generated by `multiqc` summarizing the quality trimming results |
+++