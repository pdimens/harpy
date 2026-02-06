# new
### modules
- `diagnose` now has 2 subcommands:
  - `stall`: same as previous `diagnose` behavior, where it runs snakemake with `--dry-run --debug-dag`
  - `rule`: attempt to directly run the failing rule of a workflow as identified in the snakemake log, will attempt to run snakemake to generate missing inputs if necessary
- `phase bam` added for more fine-tuned and configurable alignment phasing
- `report` added to render new ipython reports as a MySTmd website
### options
- `resume` has new `--direct` option to call Snakemake directly without harpy intervention
- hidden common option `--clean` with the options `w`, `s`, and/or `l`, to remove the `workflow/`, `.snakemake/`, and/or `logs/` directories in the output
  - this option is **hidden** because it's meant more for advanced users or development
  - options provided as sequential letters (e.g. `ws`, `sl`, `lw`, etc.)
### misc
- progress bar has a new column to show a count of the active jobs!
- time elapsed column in progress bar pauses when there are no active jobs for that rule (better reflecting the actual time elapsed)
- output log of checks and validations printed to console so Harpy is transparent about any observed delays before kicking off Snakemake
  - disabled when `--quiet` > 0
- added "workflow setup complete" text when using `--setup`
- added "all stuff is there" equivalent text when snakemake reports there is nothing to do


# changes
## breaking
### deprecations
  - `convert` (replaced by `Djinn` software)
  - `downsample` (replaced by `Djinn` software)
  - `simulate` (replaced by `Mimick` and `VISOR-HACks` softwares)
### renamed
- `demultiplex` is now `preprocess` to better reflect what the commands do
- `diagnose` is now `diagnose stall` to accomodate distinction from new `diagnose rule`
- `phase` has been renamed `phase snp` to accommodate a disctinction from the new `phase bam` workflow
- `--setup-only` replaced with the more succinct `--setup`
- the `workflow.yaml` files now all have a standard/consistent format with three main sections whose names are capitalized (whereas all the rest are lowercase):
  - `Workflow`: with common information (name, linkedread info, report skip/contigs, harpy-specific snakemake things)
  - `Parameters`: the run configurations resulting from command-line arguments/options
  - `Inputs`: the input files
  - this means previous `workflow.yaml` files are incompatible with this and future versions
- all `workflow.yaml` keys use hyphens instead of underscores to reduce keystrokes
  - e.g. `min_len` => `min-len`
- minimum mapping quality has been consistenlty named `min-map-quality` in all `workflow.yaml` files 
### function
- `sv naibr` no longer phases input alignment files, use the new `phase bam` module for that
### reports
Reports have been completely rewritten (for the third time), moving away from R/Quarto to Python/Jupyter. This change was necessary to achieve specific quality-of-life improvements that were not possible with the current setup:
- report generation is significantly faster because Quarto isn't trying to render each as a standalone HTML file
- IPython notebooks store the results within themselves, which can be accessed whenever via JupyterLab/VScode/etc
- all harpy-generated (non-MultiQC) reports can be bound and built into a singular report website using `harpy report` with navigation between pages, searching, etc.
- plots are now generated with Altair (Vega/Vegalite) plotting library, which is very fast, interactive, and reactive, with a more permissible license than Highcharts, which use a restrictive software license
### internal
- significant rewrite of the `Workflow` class and how it expects workflow, parameter, and input delcarations
- 4 SV reports consolidated into 1
- printing functions consolidated into `HarpyPrint` class

## non-breaking
- statusbar when downloading/installing workflow dependencies now lists the environment being downloaded/installed instead of saying "working..."
- progress bar does not disappear with default `--quiet` setting
- minor progress bar tweaks 
- workflow start and end time moved to row in each of the output tables rather than appear inline
- no more custom snakemake logfile handling (backported into `v3.2`)
  - `harpy view log` now points to the `.snakemake/log` folder, but is otherwise the same
- snakefiles have harpy version at the top
- workflow errors now include harpy version

### internal
- swapped order of validations/checks
  - CLI input validations are for fast basic checks (e.g. naming conventions, presence/absence)
  - harpy validations are for more involved checks (formatting, consistency, inter-parameter checks)
- Snakemake process monitoring saw a significant rewrite
  - the new internals are easier to develop and should hopefully have more consistent exiting behavior
- `resume` logic reorganized updated to match new workflow configuration design
- the `--workflow-profile` part of the snakemake command (when using hpc) has been moved to `config.yaml` to further reduce the length of the snakemake call
- grouped and single-sample leviathan variant calling now use a single consolidated snakefile
- grouped and single-sample naibr variant calling now use a single consolidated snakefile
- workflow info printing handled differently to be more flexible

# fixes
- removed redundant validations between CLI checks and harpy checks
- wording improvements for errors and doc text
- [hopefully] no more double-printing of Snakemake errors
