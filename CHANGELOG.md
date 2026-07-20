# New
## QC
- replace multiqc report with native harpy report

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