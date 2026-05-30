# Fixes
- `harpy report` error when not used in a git-enable repository
- embedded image viewer has cleaned up HTML
- no more double-printing worklow name when using `harpy resume`

# Features
## Impute
- the `--buffer` option includes fractional scaling when <= 1
  - e.g. `--strategy window:100000 -b 0.1` equivalent to `--strategy window:100000 -b 10000`
  - e.g. `--strategy contig1:1-100000 -b 0.1` equivalent to `--strategy contig1:1-100000 -b 10000`
  - only applies to window and region strategies

## Preprocess GIH
- adds `adapters.fasta` output to be used as input into `harpy qc`

# Changes
## Reports
- tables have "Export CSV" replaced with "⤓ Download CSV" to be more obvious

# Breaking
None