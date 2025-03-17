# Changes in v2.0

## New
### Non-linked-read compatibility
- `--ignore-bx` added to `harpy qc` and `harpy align` to skip linked-read routines and make the workflows suitable for non linked-read data
### Modules
- `harpy diagnose` to view snakemake's `--debug-diag` output for troubleshooting
- `harpy template` to create template files (see breaking changes below)
### Options
- `harpy sv leviathan` adds `--duplicates` and `--sharing-thresholds` options
- `harpy demultiplex gen1` adds `--keep-unknown` to retain reads that failed to demultiplex in a separate file
- `harpy demultiplex gen1` adds `--qxrx` to include the `QX:Z` and `RX:Z` tags in the read headers (defaults to not doing that)
- `harpy downsample` adds `--hpc` for cluster submission
- `harpy simulate linkedreads` adds `--merge-haplotypes` as a convenience features to merge R1 reads for hap0 and hap1 (same for R2)

## Breaking Changes
### Renamed Commands
- `harpy popgroup` is now `harpy template groupings`
- `harpy imputeparams` is now `harpy template impute`
- `harpy hpc` is now `harpy template hpc-*`
  - where `*` is the type of scheduler, e.g. `harpy template hpc-slurm`
### HPC support
- `harpy template hpc-*` now outputs `hpc/system.yaml` rather than `hpc/system/config.yaml`
- the `--hpc` option for workflows now takes the configuration yaml file directly, directory input is disabled
    - a copy is made to `outdir/workflow/hpc/config.yaml` and the snakemake workflow is configured to run off of that copy
    - input yaml file is tested for proper yaml syntax
### --quiet
- `--quiet` has been reworked as a number between 0 and 2 (i.e. `0`, `1`, `2`)
  - `0` prints starting text and progress bars as is normal (default)
  - `1` prints starting text and single progress bar for `Total`
  - `2` prints nothing, like the original `--quiet` behavior
### misc
- `harpy qc -d` parameter no longer needs commas, e.g. `harpy qc -d 10 12 14 51 -a auto ...`
- updates to `click` (internal) mean you need to call the docstring up deliberately with `harpy XXXX --help`
  - empty module call no longer brings up the docstring
- direct HTCondor support is gone in `harpy template hpc-` because the snakemake plugin seems to have vanished

## Non-breaking changes
- the molecule distance parameter for `align bwa` and `align strobe` now defaults to `0` to disable distance-based deconvolution
- quarto-based html reports now use the Lumen theme again, tweaked to look like how it did when they were originally made by FlexDashboard.
  - The theming is pulled from github upon workflow execution, so the report styling should always be up to date regardless of harpy version from now on
  - favicon is back!
- `harpy --help` docstring organized into three sections now (up from two)
- colorized the boxes for docstrings
- duplicate-marking in the `align` workflows now internally sets a distance threshold for determining optical duplicates based on sequencing platform
- coverage and molecular coverage are calculated faster and more accurately
- `demultiplex gen1` now uses the purpose-built `dmox` for significantly faster demultiplexing
- different spinner for progress bars
- the BUSCO analysis in `harpy assembly` and `harpy metassembly` is now hardcoded to use orthoDB v12
- the unified `Genome/` folder has been replaced with a per-output-folder `outdir/workflow/genome/` folder that serves the same function to avoid situations where multiple fasta files with the same name are being used simultaneously across concurrent workflows.
  - the lack of unification results in extra redundancy (more disk space used =/) but with the benefit of more reliable behavior 

## Fixes
- improved error handling in quarto-based reports when there are NAs or NaNs, or too few data points to make a histogram
- threads are capped at `999` for workflows to prevent Snakemake from getting in your face about it
- html injection into reports with circos plots to make sure the tooltip isn't hidden behind the plot
- bugfixes in the `lsf` profile from `harpy hpc lsf`
- report fixes for plot sizes, axes limits, etc.
- `harpy view` uses better internal logic
  - prevents the "broken pipe" message after closing a gzipped file
  - no more [direct] reliance on `subprocess` calls
  - calls `pygmentize` directly through python API
- updated mentions of `popgroup` and `stitchparams` (legacy) to `harpy template`

## Internal
- `popgroup.py` and `imputeparams.py` merged into `template.py`
- the `click` bindings for `simulate` were moved out of `__main__.py` and into `simulate.py`
  - result is less cluttered and more intuitive `__main__.py`