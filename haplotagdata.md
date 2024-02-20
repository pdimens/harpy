---
label: Haplotag data
icon: file-binary
order: 5
---
# :icon-file-binary: Haplotag data

## Data Format
### Barcodes
While barcodes are actually combinatorial bases, in the read headers they are represented
with the format `AxxCxxBxxDxx`, where each barcode segment is denoted as `Axx` (or `Bxx`, etc.).
The capital letter denotes which preparation microwell plate the barcode segment is from (plate `A`, `B`, `C`, or `D`) 
and `xx` is a number between `00` and `96` corresponding to the well from that microplate.
So, the `A14` segment would correspond with the barcode from Plate `A`, well `14`.
A `00` barcode (e.g. `C00`) indicates a missing/invalid barcode segment, which invalidates the entire barcode.

#### barcode protocol varieties
If you think haplotagging is as simple as exactly $96^4$ unique barcodes, you would only be half-correct. The original haplotagging
protocol in Meier *et al.* is good, but the authors (and others) have been working to improve this linked-read technology to improve
things like reduce PCR duplicates, improve successful barcode sequencing and error correction, etc. As a result, a few updated variants
of the beadtags will appear, each likely with their own way of properly demultiplexing the raw sequences. Harpy will aim to incorporate
demultiplexing these beadtag variants as they become available.

#### where the barcodes go
Chromium 10X linked-reads have a particular format where the barcode is the leading 16 bases 
of the read. However, haplotagging data **does not use that format**, nor do the tools 
implemented in Harpy work correctly with it. Once demultiplexed, haplotagging sequences should look 
like regular FASTQ files of inserts and the barcode is stored in a `BX:Z:AxxCxxBxxDxx` tag 
in the read header. Again, **do not include the barcode in the sequence**.

### Read headers
Like mentioned, the haplotag barcode is expected to be stored in the `BX:Z:` tag in the 
read header. This information is retained through the various Harpy
steps. An example read header could look like:
``` example valid read header
@A00814:267:HTMH3DRXX:2:1101:4580:1000  RX:Z:GAAACGACCAACA+CGAACACGTTAGC    QX:Z:F,FFFFFFFFFFF+FF,FFF:FFFFFF   BX:Z:A02C01B11D46
```
Notably, only the sequence ID (`@...`) and `BX:Z:` tag are actually required. In the example 
above, there are additional tags (`RX:Z:` and `QX:Z:`) which arent used by Harpy, but they 
conform to the [SAM comment spec (section 1.5)](https://samtools.github.io/hts-specs/SAMv1.pdf) 
of `TAG:TYPE:VALUE`. The takeaway is that the `BX:Z:` tag can be anywhere in the read header 
after the sequence ID as long as any tags after it conform to the SAM spec `TAG:TYPE:VALUE` (see note). 
Read comments that aren't following the `TAG:TYPE:VALUE` SAM spec may cause downstream errors.  

!!!warning A caveat
The Leviathan structural variant caller expects the `BX:Z:` tag at the end of the alignment 
record, so if you intend on using that variant caller, you will need to make sure the `BX:Z:`
tag is the last one in the _sequence alignment_ (BAM file). If you use Harpy to align the 
sequences, then it will make sure the `BX:Z:` tag is moved to the end of the alignment.
!!!

### Read length
Reads must be at least 30 base pairs in length for alignment. The `qc` module removes reads <50bp.

### Compression
Harpy generally doesn't require the input sequences to be in gzipped/bgzipped format, but it's good practice to compress your reads anyway.
Compressed files are expected to end with the extension `.gz`.

### Naming conventions
Unfortunately, there are many different ways of naming FASTQ files, which makes it 
difficult to accomodate every wacky iteration currently in circulation.
While Harpy tries its best to be flexible, there are limitations. 
To that end, for the `demultiplex`, `qc`, and `align` modules, the 
most common FASTQ naming styles are supported:
- **sample names**: Alphanumeric and `.`, `-`, `_`
    - you can mix and match special characters, but that's bad practice and not recommended
    - examples: `Sample.001`, `Sample_001_year4`, `Sample-001_population1.year2` <- not recommended
- **forward/reverse**: `_F`, `.F`, `_R1`, `.R1`, `_R1_001`, `.R1_001`, *etc.*
    - note that this **does not include** `.1` or `_1` conventions for forward/reverse
- **fastq extension**: `.fastq`, `.FASTQ`, `.fq`, `.FQ`
- **gzipped**: supported and recommended
- **not gzipped**: supported

You can also mix and match different formats and styles within a given directory, although again, **this isn't recommended**.
As a good rule of thumb for any computational work, you should be deliberate and consistent in how you name things.

## Barcode thresholds
By the nature of linked read technologies, there will (almost always) be more DNA fragments than unique barcodes for them. As a result,
it's common for barcodes to reappear in sequences. Rather than incorrectly assume that all sequences/alignments with the same barcode
originated from the same orignal DNA molecule, linked-read aware programs include a threshold parameter to determine a "cutoff" distance
between alignments with the same barcode. This parameter can be interpreted as "if a barcode appears more than `X` base pairs away from the
same barcode (on the same contig), then we'll consider them as originating from different molecules." If this threshold is lower, then
you are being more strict and indicating that alignments sharing barcodes must be closer together to be considered originating from the same
DNA molecule. Conversely, a higher threshold indicates you are being more lax and indicating barcodes can be further away from each other
and still be considered originating from the same DNA molecule. A threshold of 50kb-150kb is considered a decent balance, but you should choose
larger/smaller values if you have evidence to support them. 

![Molecule origin is determined by the distance between alignments with the same barcode relative to the specified threshold](/static/bc_threshold.png)

| Alignment distance     |    Inferred origin  |
|:-----------------------|:--------------------|
| less than threshold    |     same molecule   |
| greater than threshold | different molecules |