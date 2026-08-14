# New
## SNP
- `-i`/`--keep-invariant` to retain monomorphic sites

## Improvements
- variant sorting of temporary region-chunks in `snp` workflows now outputs uncompressed BCF, which should see a bit of speedup
- `preprocess gih` stagger code should be a little leaner with fewer allocations and no BAM->SAM conversion

## Fixes
- report tables might respect dark mode better. might not. I'm not a web dev :shrug: