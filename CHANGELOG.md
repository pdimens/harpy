# Fixes
## QC
- adapter trimming was being incorrectly skipped in all cases
- set max threads for fastp to 4
## Align
- both `strobealign` and `bwa` workflows have the rules shifted a bit for better memory management
- `harpy-utils optical-distance` now properly falls back to `100`
- `harpy-utils molecule-coverage` is much less RAM hungry
- `harpy-utils bx-stats-sam` is much less RAM hungry
