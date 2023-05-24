---
label: Haplotag data format
icon: file-binary
order: 97
---
# :icon-file-binary: Haplotag data format

## Barcodes
While barcodes are actually combinatorial bases, in the read headers they are represented
with the format `AxxCxxBxxDxx`, where each barcode segment is denoted as `Axx` (or `Bxx`, etc.).
The capital letter denotes which preparation microwell plate the barcode segment is from (plate `A`, `B`, `C`, or `D`) 
and `xx` is a number between `00` and `96` corresponding to the well from that microplate.
So, the `A14` segment would correspond with the barcode from Plate `A`, well `14`.
A `00` barcode (e.g. `C00`) indicates a missing/invalid barcode segment, which invalidates the entire barcode.

### where the barcodes go
Chromium 10X linked-reads have a particular format where the barcode is the leading 16 bases 
of the read. However, haplotagging data **does not use that format**, nor do the tools 
implemented in Harpy work correctly with it. Simply put, haplotagging sequences should look like regular FASTQ files of inserts and the barcode is stored in a `BX:Z:` tag in the read header. Again, **do not include the barcode in the sequence**.

## Read headers
Like mentioned immediately above this sentence, the haplotag barcode is expected to be stored
in the `BX:Z:` tag in the read header. This information is retained through the various Harpy
steps. An example read header could look like:
``` example valid read header
@A00814:267:HTMH3DRXX:2:1101:4580:1000 BX:Z:A02C01B11D46        RX:Z:GAAACGACCAACA+CGAACACGTTAGC   QX:Z:F,FFFFFFFFFFF+FF,FFF:FFFFFF
```
Notably, only the sequence ID (`@...`) and `BX:Z:` tag are actually required. In the example 
above, there are additional tags (`RX:Z:` and `QX:Z:`) which arent used by Harpy, but they 
conform to the [SAM comment spec (section 1.5)](https://samtools.github.io/hts-specs/SAMv1.pdf) 
of `TAG:TYPE:VALUE`. The takeaway is that the `BX:Z:` tag can be anywhere in the read header 
after the sequence ID as long as any tags after it conform to the SAM spec `TAG:TYPE:VALUE`. Read comments that aren't following the `TAG:TYPE:VALUE` SAM spec may cause downstream errors (definitely prevents LEVIATHAN from running).  

## Read length
Reads must be at least 30 base pairs in length for alignment. The `trim` module removes reads <50bp.

## Compression
Harpy expects the input sequences to be in gzipped/bgzipped format, therefore the file names should end in `.gz`. This is a hard requirement and it's good practice to compress your reads
anyway. 

## Naming conventions
Input sequences will have the format `{samplename}.{extension}`. To make sure there are no hiccups, this section details valid naming conventions for both the `{samplename}` and `{extention}` parts of file names intended as input for Harpy.

### sample names
Sample names must not have any unusual special characters or dot (`.`) separators in them.
- 🔵 good: `sample_01_pop1`
    - example: `sample_01_pop1.F.fastq.gz`
- 🔵 good: `sample-01`
    - example: `sample-01.R.fastq.gz`
- 🚫 **bad**: `sample.01.pop1`
    - example: `sample_01.pop1.F.fastq.gz`

### file extensions
There are a handful of "accepted" naming schemes for fastq file extensions, but Harpy only 
accepts a limited number of them, shown below. The fastq files **must be consistent** with regards to the extensions and read-pair naming styles.
That is, you must only use `.fastq.gz` or only use `.fq.gz` for all files, and the same for `.
R1.`/`.R2.` or `_R1.`/`_R2.` (adhere to a single row in the tables below).
Notice that the read pair part differs from the accepted fastq names for read trimming.

==- :icon-check-circle: acceptable fastq names for **trimming**

| forward-reverse notation | extension  | example forward          | example reverse         |
|:-------------------------|:-----------|:-------------------------|:------------------------|
| `.F` / `.R`                | `fastq.gz` | ` samplename.F.fastq.gz` | `samplename.R.fastq.gz` |
| `.F` / `.R`                | `fq.gz`    | `samplename.F.fq.gz`     | `samplename.R.fq.gz`    |
| `.1` / `.2`                | `fastq.gz` | `samplename.2.fastq.gz`  | `samplename.2.fastq.gz` |
| `.1` / `.2`                | `fq.gz`    | `samplename.1.fq.gz`     | `samplename.2.fq.gz`    |

==- :icon-check-circle: acceptable fastq names for **aligning**
Notice that the forward/reverse notation is slightly different for reads expected to be used for
alignment. This is deliberate to make sure trimmed reads have a different convention from raw reads and avoid confusion (i.e. "has this file already been trimmed??").

| forward-reverse notation | extension  | example forward           | example reverse          |
|:-------------------------|:-----------|:--------------------------|:-------------------------|
| `.R1` / `.R2`            | `fastq.gz` | ` samplename.R1.fastq.gz` | `samplename.R2.fastq.gz` |
| `.R1` / `.R2`            | `fq.gz`    | `samplename.R1.fq.gz`     | `samplename.R2.fq.gz`    |
| `_R1` / `_R2`            | `fastq.gz` | `samplename_R1.fastq.gz`  | `samplename_R2.fastq.gz` |
| `_R1` / `_R2`            | `fq.gz`    | `samplename_R1.fq.gz`     | `samplename_R2.fq.gz`    |
===
