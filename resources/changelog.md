# new
- `diagnose` now has 2 subcommands:
  - `stall`: same as previous `diagnose` behavior, where it runs snakemake with `--dry-run --debug-dag`
  - `rule`: attempt to directly run the failing rule of a workflow as identified in the snakemake log, will attempt to run snakemake to generate missing inputs if necessary
- `harpy resume` has new `--direct` option to call Snakemake directly without harpy intervention
- hidden common option `--clean` with the options `w`, `s`, and/or `l`, to remove the `workflow/`, `.snakemake/`, and/or `logs/` directories in the output
  - this option is **hidden** because it's meant more for advanced users or development
  - options provided as sequential letters (e.g. `ws`, `sl`, `lw`, etc.)

# deprecations
- harpy convert (replaced by Djinn)
- harpy downsample (replaced by Djinn)
- harpy simulate (replaced by Mimick and VISOR-HACks)

# changes
- added "workflow setup complete" text when using `--setup-only`

# fixes


# In progress but not ready
- replace manually hacked progressbars with executor plugin