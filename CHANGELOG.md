# Fixes
## QC
- adapter trimming was being incorrectly skipped in all cases
- set max threads for fastp jobs to 4
## Align
- both `strobealign` and `bwa` workflows have the rules shifted for better memory management
    - specifically, the pipes between align->fixmates->markdups have been broken for speed reasons (it was very slow)
    - the new rules create temporary (uncompressed BAM) files at the choke-point steps (i.e. collate, sort, markdups) to make better use of time and computational resources
    - `samtools sort` is now given more threads and RAM per thread, with the memory per thread decreasing on failed attempts
        - 3 total attempts. Initially 3GB RAM per thread (x4 threads), drops by half each attempt (e.g. 12GB, 6GB, 3GB total)
- `samtools stats` properly ignores duplicates on processed alignments
- `harpy-utils optical-distance` now properly falls back to `100`
- `harpy-utils molecule-coverage` is much less RAM hungry
- `harpy-utils bx-stats-sam` is much less RAM hungry
## Assembly/Metassembly
- BUSCO now scrubs the orthodb version from the expected output `.txt` file to prevent orthodb version updates breaking the workflow