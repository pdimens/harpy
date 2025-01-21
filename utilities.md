---
label: Utilities
icon: terminal
order: 2
---

# :icon-terminal: Utilities
Harpy is the sum of its parts and some of those parts are [stand-alone scripts](https://github.com/pdimens/harpy/tree/main/harpy/bin)
used by the workflows that are accessible from within the Harpy conda environment.
This page serves to document those scripts, since using them outside of a workflow
might be useful too. You can call up the docstring for any one of these utilities
by calling the program without any arguments.

### assign_mi.py
```bash
assign_mi.py -c cutoff input.bam > output.bam
```
Assign an `MI:i` (Molecular Identifier) tag to each barcoded
record based on a molecular distance cutoff. Input file **must be coordinate sorted**.
This is similar to [deconvolve_alignments.py](#deconvolve_alignmentspy), except it does not record the deconvolution in the `BX` tag.
- unmapped records are discarded
- records without a `BX:Z` tag or with an invalid barcode (`00` as one of its segments) are presevered but are not assigned an `MI:i` tag

### bx_stats.py
```bash
bx_stats.py input.bam > output.gz
```
Calculates various linked-read molecule metrics from the (coordinate-sorted) input alignment file.
Metrics include (per molecule): 
- number of reads
- position start
- position end
- length of molecule inferred from alignments
- total aligned basepairs
- total length of inferred inserts
- molecule coverage (%) based on aligned bases
- molecule coverage (%) based on total inferred insert length

### bx_to_end.py
```bash
bx_to_end.py input.[fq|bam] > output.[fq.gz|bam]
```
Parses the records of a FASTQ or BAM file and moves the `BX:Z` tag, if present, to
the end of the record, which makes the data play nice with LRez/LEVIATHAN. During 
alignment, Harpy will automatically move the `BX:Z` tag to the end of the alignment
record, so that will **not** require manual intervention.

### check_bam.py
```bash
check_bam.py input.bam > output.txt
```
Parses an aligment file to check:
- if the sample name matches the `RG` tag
- whether `BX:Z` is the last tag in the record
- the counts of: 
    - total alignments
    - alignments with an `MI:i` tag
    - alignments without `BX:Z` tag
    - incorrect `BX:Z` tag

### check_fastq.py
```bash
check_bam.py input.bam > output.txt
```
Parses a FASTQ file to check if any sequences don't conform to the SAM spec,
whether BX:Z: is the last tag in the record, and the counts of:
- total reads
- reads without `BX:Z` tag
- reads with incorrect `BX:Z` tag

### concatenate_bam.py
```bash
concatenate_bam.py [--bx] file_1.bam file_2.bam...file_N.bam > output.bam
# or #
concatenate_bam.py [--bx] -b bam_files.txt > output.bam
```
Concatenate records from haplotagged SAM/BAM files while making sure `MI` tags  remain unique for every sample.
This is a means of accomplishing the same as `samtools cat`, except all `MI` tags are updated
so individuals don't have overlapping `MI` tags (which would mess up all the linked-read data). You can either provide
all the files you want to concatenate, or a single file featuring filenames with the `-b` option. Use the `--bx` option
to also rewrite `BX` tags such that they are unique for every individual too, although take note that there can only be
$96^4$ (84,934,656) unique haplotag barcodes and it will raise an error if that number is exceeded.
 
### count_bx.py
```bash
count_bx.py input.fastq > output.txt
```
Parses a FASTQ file to count:
- total sequences
- total number of `BX` tags
- number of valid haplotagging `BX` tags
- number of invalid `BX` tags
- number of invalid `BX` tag segments (i.e. `A00`, `C00`, `B00`, `D00`).

### deconvolve_alignments.py
```bash
deconvolve_alignments.py -c cutoff input.bam > output.bam
```
Deconvolve BX-tagged barcodes and assign an `MI` (Molecular Identifier) tag to each barcoded record based on a molecular distance cutoff.
Input file **must be coordinate sorted**. This is similar to [assign_mi.py](#assign_mipy), except it will also deconvolve the `BX` tag by
hyphenating it with an integer (e.g. `A01C25B31D92-2`).
- unmapped records are discarded
- records without a `BX` tag or with an invalid barcode (`00` as one of its segments) are presevered but are not assigned an `MI` tag

### depth_windows.py
```bash
samtools depth -a file.bam | depth_windows.py windowsize > output.txt
```
Reads the output of `samtools depth -a` from stdin and calculates means within windows of a given `windowsize`.

### haplotag_acbd.py
```bash
haplotag_acbd.py output_directory
```
Generates the `BC_{ABCD}.txt` files necessary to demultiplex Gen I haplotag barcodes into the specified `output_directory`.

### infer_sv.py
```bash
infer_sv.py file.bedpe [-f fail.bedpe] > outfile.bedpe
```
Create column in NAIBR bedpe output inferring the SV type from the orientation. Removes variants with FAIL flags
and you can use the optional `-f` (`--fail`) argument to output FAIL variants to a separate file.

### inline_to_haplotag.py
```bash
inline_to_haplotag.py -f <forward.fq.gz> -r <reverse.fq.gz> -b <barcodes.txt> -p <prefix> > barcodes.conversion.txt
```
Converts inline nucleotide barcodes in reads to haplotag linked reads with barcodes in `BX:Z` and `OX:Z` header tags.

### leviathan_bx_shim.py
```bash
leviathan_bx_shim.py input.bam > output.bam
```
Uses the `MI` tags in `input.bam` as a point of reference to deconvolve the `BX` tags by rewriting them as unique,
non-hyphenated ACBD tags. This "shim" script is necessary to preprocess a BAM file prior to variant calling in
LEVIATHAN because the software relies exclusively on `BX` tags and does not perform any internal deconvolution. 
- requires `input.bam` to have `MI` tags for alignment records
  - all `harpy align` workflows provide `MI` tags, but these may not be present if the BAM file was generated by other means.

### make_windows.py
```bash
make_windows.py -w <window.size> -m <0,1> input.fasta[.fai] > output.bed
```
Create a BED file of fixed intervals (`-w`, --`window`) from a FASTA or fai file (the kind generated with `samtools faidx`).
Nearly identical to `bedtools makewindows`, except the intervals are nonoverlapping. The `-m` (`--mode`) option specified
whether indexing starts at `0` or `1`.

### molecule_coverage.py
```bash
molecule_coverage.py -f genome.fasta.fai statsfile > output.cov
```
Using the statsfile generated by `bx_stats.py` from Harpy, will calculate "molecular coverage" across the genome.
Molecular coverage is the "effective" alignment coverage if you treat a molecule inferred from linked-read data as
one contiguous alignment, even though the reads that make up that molecule don't cover its entire length. Requires a
FASTA fai index (the kind created with `samtools faidx`) to know the actual sizes of the contigs.

### parse_phaseblocks.py
```bash
parse_phaseblocks.py input > output.txt
```
Parse a phase block file from HapCut2 to pull out summary information

### rename_bam
```bash
rename_bam.py [-d] new_name input.bam
```
Rename a sam/bam file and modify the `@RG` tag of the alignment file to reflect the change for both `ID` and `SM`.
This process creates a new file `new_name.bam` and you may use `-d` to delete the original file. Requires `samtools`.

### separate_singletons
```bash
separate_singletons -t threads -b barcode_tag -s singletons.bam input.bam > output.bam
```
Isolate singleton and non-singleton linked-read BAM records into separate files. Singletons
refers to barcodes that have only one unpaired or paired read, meaning the barcode doesn't
actually link and reads togeher.

### separate_validbx
```bash
separate_validbx invalid.bam input.bam > valid.bam
```
Split a BAM file with `BX` tags by tag validity
- `stdout`: alignments with valid ACBD barcodes
- **first argument**: alignments invalid ACBD barcodes
