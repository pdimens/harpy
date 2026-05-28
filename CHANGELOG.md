# Fixes
- `harpy report` error when not used in a git-enable repository

# Features
## Impute
- the `--buffer` option includes fractional scaling when <= 1
  - e.g. `--strategy window:100000 -b 0.1` equivalent to `--strategy window:100000 -b 10000`
  - e.g. `--strategy contig1:1-100000 -b 0.1` equivalent to `--strategy contig1:1-100000 -b 10000`
  - only applies to window and region strategies