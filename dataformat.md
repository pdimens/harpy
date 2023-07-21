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

### barcode protocol varieties
If you think haplotagging is as simple as exactly $96^4$ unique barcodes, you would only be half-correct. The original haplotagging
protocol in Meier *et al.* is good, but the authors (and others) have been working to improve this linked-read technology to improve
things like reduce PCR duplicates, improve successful barcode sequencing and error correction, etc. As a result, a few updated variants
of the beadtags will appear, each likely with their own way of properly demultiplexing the raw sequences. Harpy will aim to incorporate
demultiplexing these beadtag variants as they become available.

### where the barcodes go
Chromium 10X linked-reads have a particular format where the barcode is the leading 16 bases 
of the read. However, haplotagging data **does not use that format**, nor do the tools 
implemented in Harpy work correctly with it. Simply put, haplotagging sequences should look like regular FASTQ files of inserts and the barcode is stored in a `BX:Z:AxxCxxBxxDxx` tag in the read header. Again, **do not include the barcode in the sequence**.

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
Harpy generally doesn't require the input sequences to be in gzipped/bgzipped format, but it's good practice to compress your reads anyway.
Compressed files are expected to end with the extension `.gz`. There is a hard requirement for fastq files to be gzipped when [using ema](#file-naming-for-align-ema) to align reads. 

## Naming conventions
Unfortunately, there are many different ways of naming FASTQ files, which makes it difficult to accomodate every wacky iteration currently in circulation.
While Harpy tries its best to be flexible, there are limitations, specifically with how Snakemake has difficulty with the `harpy align --method ema` module.
To that end, for the `demultiplex`, `trim`, and `align --method bwa` methods, the most common FASTQ naming styles are supported:
- **sample names**: Alphanumeric and `.`, `-`, `_`
    - you can mix and match special characters, but that's bad practice and not recommended
    - examples: `Sample.001`, `Sample_001_year4`, `Sample-001_population1.year2` <- not recommended
- **forward/reverse**: `_F`, `.F`, `_R1`, `.R1`, `_R1_001`, `.R1_001`, *etc.*
- **fastq extension**: `.fastq`, `.FASTQ`, `.fq`, `.FQ`
- **gzipped**: supported
- **not gzipped**: supported (but not recommended)

You can also mix and match different formats and styles within a given directory, although again, this isn't recommended
As a good rule of thumb for any computational work, you should be deliberate and consistent in how you name things.

### file naming for align-ema
As mentioned, the `align` workflow using `ema` is a bit of a problem with this. Not for lack of trying, but the logic that keeps the other modules
flexible with naming doesn't seem to work with align-ema, at least not when Snakemake is concerned. Therefore, if you want to align you reads with `ema`,
then you need to make sure your reads follow a **specific** and **consistent** naming convention. We would happily welcome anyone who is interested in fixing
that module's logic to be as flexible as the other ones 🙏. 

#### sample names
Follow best-practices-- if you need to use special characters (among `.`, `_`, and `-`), then pick one and stick with it.
When in doubt, use underscores (`_`).
- examples: `Sample001`, `Sample_001`, `Sample001_population002_year3` 

#### file extensions 
👉 TL;DR: adhere to a single row in the table below 👈

The fastq files **must be consistent** with regards to the extensions and read-pair naming styles.
That is, you must only use `.fastq.gz` or only use `.fq.gz` for all files, and the same for `.
R1.`/`.R2.` or `_R1.`/`_R2.`. Note that these are case-sensitive and the F/R1 designations **do not** have `_001` after them, like
is sometimes seen in fastq names (*e.g.,* `_R1_001.fastq.gz`).

#### :icon-check-circle: acceptable fastq names for **aligning with ema**
| forward-reverse notation | extension  | example forward           | example reverse          |
|:-------------------------|:-----------|:--------------------------|:-------------------------|
| `.R1` / `.R2`            | `.fastq.gz` | ` samplename.R1.fastq.gz` | `samplename.R2.fastq.gz` |
| `.R1` / `.R2`            | `.fq.gz`    | `samplename.R1.fq.gz`     | `samplename.R2.fq.gz`    |
| `_R1` / `_R2`            | `.fastq.gz` | `samplename_R1.fastq.gz`  | `samplename_R2.fastq.gz` |
| `_R1` / `_R2`            | `.fq.gz`    | `samplename_R1.fq.gz`     | `samplename_R2.fq.gz`    |

