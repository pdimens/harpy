# New
## QC
- replace multiqc report with native harpy report
## Align
- BWAMEM2 has been replaced with minibwa. Long live BWA!
  - it's much faster, and takes much less time to index a reference
- coverage depth added to aggregate report for processed alignments 
- minimap2 added back in for long-read compatability
  - called with `harpy align minimap`

## Report
- harpy reports can be converted to less-nice but functional standalone HTML files
  - this feature is accessed using `harpy report static`
  - to accomodate this, `harpy report` (live report website) is now `harpy report live`

# Changes
## Align
- `-d` (molecule distance threshold) has its default restored to 50kb since this value is used exclusively for reporting and does not alter the data

## misc
- removed FASTA format validation because it can be dreadfully slow with existing tools



# Fixes
## reports
- tables now render properly in VScode/Jupyter contexts

## misc
- constrain CASAVA regex in FASTQ file validation so it doesn't trigger false positives when new CASAVA appears in unexpected places
- [internal] notebooks no longer a submodule/subdirectory of `harpy.report`
- utility `check_fastq.py` no longer employs globals, instead uses a sensible class system
- add multithreading to pre-workflow VCF and XAM file validation and parsing
- [internal-ish] the bwa, strobealign, and minimap2 workflows are nearly identical except for the reference preprocessing and alignment, so to minimize redundancy and duplication, those workflows have a single `align.smk` that imports a second snakefile `align_{aligner}.smk` that handles just the preprocessing and direct alignment for those aligners, then hands off to `align.smk` for all the downstream things (dedup, sorting, reports, etc)