# New
## QC
- replace multiqc report with native harpy report
## Align
- BWAMEM2 has been replaced with minibwa. Long live BWA!
  - it's much faster, and takes much less time to index a reference
- coverage depth added to aggregate report for processed alignments 

## Report
- harpy reports can be converted to less-nice but functional standalone HTML files
  - this feature is accessed using `harpy report static`
  - to accomodate this, `harpy report` (live report website) is now `harpy report live`

# Fixes
## reports
- tables now render properly in VScode/Jupyter contexts

## misc
- constrain CASAVA regex in FASTQ file validation so it doesn't trigger false positives when new CASAVA appears in unexpected places
- [internal] notebooks no longer a submodule/subdirectory of `harpy.report`
- utility `check_fastq.py` no longer employs globals, instead uses a sensible class system
- add multithreading to pre-workflow VCF and XAM file validation and parsing