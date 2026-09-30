# New
- `harpy view envs` prints versions too, where applicable
- CLI validations for `--hpc`/`-H` to check if the plugins necessary for the configuration are installed
- `harpy preprocess meier2021` now uses `dmox` v0.3.1
  - introduces the `--use-stitch-base` logic via `--stitch`
  - introduces the `--allow-complementary-stitch` logic via `--stitch-comp`
- `harpy align arachne` adds Arachne linked-read aware aligner (successor to lariat)
- tab-completion for `harpy` in bash, zsh, and fish
  - the completion scripts are generated during the conda/pixi build, installed where each shell looks for them, and enabled automatically when the environment is activated (conda, `pixi shell`), no user setup required
  - file-like arguments (FASTQ, FASTA, BAM, VCF, HPC profiles, etc.) defer to the shell's own file completion
  - adds the hidden `harpy completion` command that prints the script for a given shell (used by the build)

# Changes
- much faster CLI startup: `harpy --version` and `harpy --help` go from ~2s to ~0.2s
  - the subcommands of `harpy` and `harpy-utils` are now imported only when they are run
  - heavy imports (e.g. `nbconvert`, `pysam`, `yaml`) moved out of module scope into the functions that need them
- some workflows with large temporary files (like `align`) have jobs grouped to prioritize running steps that would:
  1. remove the temporary file sooner
  2. if running on an HPC, reduces the need to copy temp files between nodes
- the `dmox` update makes `workflow.yaml` files from previous harpy versions invalid
  - not technically invalid, but the key absence/mismatch will default to all optional features turned off, which may be unintended
- fastq validation is now limited to 100 records, which should see a significant validation speedup

# Fixes
- `harpy view envs`: simpler logic and print diagnostic text if empty
- more robust snakemake error printing (again)
  - this time it's a parse-and-gather approach that uses an internal class
  - strengthened outputting snakemake missing and syntax errors
- `harpy resume` no longer overwrites the harpy version of `workflow.yaml`
- error printing when using `--container` correctly displays full apptainer-prefixed shell call


# Internal
- simplified summaries logic
- harpy-utils: replace Golang `regexp` with `coregex` for speed/efficiency
- progressbar logic moved out of `HarpyPrint` class in `printing.py` and into `progress.py`
  - removed the redundant `rich.Live` wrap to the progress bars-- progress bars should be a little snappier
- swapped pandas for polars (speed!)


# Added but not exposed
- `harpy snp deepvariant` - call SNPs in high-depth samples using Google's DeepVariant. AI AI AI!!!
  - this option must use the docker image of the software because of the way DeepVariant is packaged, so `--container` isn't exposed
- this is just internal scaffolding for now. Can be made public with sufficient interest