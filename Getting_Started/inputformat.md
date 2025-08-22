---
label: Input Format
icon: file-binary
order: 100
---
# :icon-file-binary: Input Format

## Read length
Reads must be at least 30 base pairs in length for alignment. By default, the [!badge corners="pill" text="qc"](/Workflows/qc.md) module removes reads <30bp.

## Compression
Harpy generally doesn't require the input sequences to be in gzipped/bgzipped format, but it's good practice to compress your reads anyway.
Compressed files are expected to end with the extension [!badge variant="success" text=".gz"].

## Naming conventions
Unfortunately, there are many different ways of naming FASTQ files, which makes it 
difficult to accomodate every wacky iteration currently in circulation.
While Harpy tries its best to be flexible, there are limitations. 
To that end, for the [!badge corners="pill" text="deumultiplex"](/Workflows/demultiplex.md), [!badge corners="pill" text="qc"](/Workflows/qc.md), and [!badge corners="pill" text="align"](/Workflows/Align/bwa.md) modules, the 
most common FASTQ naming styles are supported:
- **sample names**: [!badge variant="success" text="a-z"] [!badge variant="success" text="0-9"] [!badge variant="success" text="."] [!badge variant="success" text="_"] [!badge variant="success" text="-"] [!badge variant="secondary" text="case insensitive"]
    - you can mix and match special characters, but that's bad practice and not recommended
    - examples: `Sample.001`, `Sample_001_year4`, `Sample-001_population1.year2` <- not recommended
- **forward**: [!badge variant="success" text="_F"] [!badge variant="success" text=".F"] [!badge variant="success" text="_1"] [!badge variant="success" text=".1"] [!badge variant="success" text="_R1_001"] [!badge variant="success" text=".R1_001"] [!badge variant="success" text="_R1"] [!badge variant="success" text=".R1"] 
- **reverse**: [!badge variant="success" text="_R"] [!badge variant="success" text=".R"] [!badge variant="success" text="_2"] [!badge variant="success" text=".2"] [!badge variant="success" text="_R2_001"] [!badge variant="success" text=".R2_001"] [!badge variant="success" text="_R2"] [!badge variant="success" text=".R2"] 
- **fastq extension**: [!badge variant="success" text=".fq"] [!badge variant="success" text=".fastq"] [!badge variant="secondary" text="case insensitive"]
- **gzipped**: [!badge variant="info" icon=":thumbsup:" text="supported"] [!badge variant="info" icon=":heart:" text="recommended"]
- **not gzipped**: [!badge variant="info" icon=":thumbsup:" text="supported"]

You can also mix and match different formats and styles within a given directory, although again, **this isn't recommended**.
As a good rule of thumb for any computational work, you should be deliberate and consistent in how you name things.

