## new
- progress bar has a new column to show a count of the active jobs!
- time elapsed column in progress bar pauses when there are no active jobs for that rule (better reflecting the actual time elapsed)
- `diagnose` now has 2 subcommands:
  - `stall`: same as previous `diagnose` behavior, where it runs snakemake with `--dry-run --debug-dag`
  - `rule`: attempt to directly run the failing rule of a workflow as identified in the snakemake log, will attempt to run snakemake to generate missing inputs if necessary
- `harpy resume` has new `--direct` option to call Snakemake directly without harpy intervention
- hidden common option `--clean` with the options `w`, `s`, and/or `l`, to remove the `workflow/`, `.snakemake/`, and/or `logs/` directories in the output
  - this option is **hidden** because it's meant more for advanced users or development
  - options provided as sequential letters (e.g. `ws`, `sl`, `lw`, etc.)
- output log of checks and validations printed to console so Harpy is transparent about the delays before kicking off Snakemake
  - disabled when `--quiet` > 0
- added "workflow setup complete" text when using `--setup`

## changes
### breaking
- deprecations
  - harpy convert (replaced by Djinn)
  - harpy downsample (replaced by Djinn)
  - harpy simulate (replaced by Mimick and VISOR-HACks)
- `--setup-only` replaced with the more succinct `--setup`

### non-breaking
- swapped order of validations/checks
  - CLI input validations are for fast basic checks (e.g. naming conventions, presence/absence)
  - harpy validations are for more involved checks (formatting, consistency, inter-parameter checks)
- [internal] Snakemake process monitoring saw a significant rewrite
  - the new internals are easier to develop and should hopefully have more consistent exiting behavior
- [internal] `harpy resume` logic reorganized a bit

# fixes
- removed redundant validations between CLI checks and harpy checks
- wording improvements for errors and doc text
- [hopefully] no more double-printing of Snakemake errors
