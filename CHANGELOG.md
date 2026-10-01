# New
- `harpy view envs` prints versions too, where applicable
- CLI validations for `--hpc`/`-H` to check if the plugins necessary for the configuration are installed
- `harpy preprocess meier2021` now uses `dmox` v0.3.1
  - introduces the `--use-stitch-base` logic via `--stitch`
  - introduces the `--allow-complementary-stitch` logic via `--stitch-comp`
- `harpy align arachne` adds Arachne linked-read aware aligner (successor to lariat)
- tab-completion for `harpy`, `harpy-utils`, and `hv` in bash, zsh, and fish
  - the completion scripts are generated during the conda/pixi build, installed where each shell looks for them, and enabled automatically when the environment is activated (conda, `pixi shell`), no user setup required
  - file-like arguments (FASTQ, FASTA, BAM, VCF, HPC profiles, etc.) defer to the shell's own file completion
  - adds the hidden `harpy completion <shell> [program]` command that prints the script for a given shell and program (used by the build)

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
- `bx-stats-sam` correctly names the column `fragments`, was formerly `reads`

# Fixes
- `harpy view envs`: simpler logic and print diagnostic text if empty
- more robust snakemake error printing (again)
  - this time it's a parse-and-gather approach that uses an internal class
  - strengthened outputting snakemake missing and syntax errors
  - workflow errors when running with an HPC scheduler (`--hpc`) now show what can be found in snakemake's output: the failing rule, snakemake's message, the scheduler's job ID, and the contents of the job's log files (including the scheduler's own log for slurm, lsf, and googlebatch)
    - log files snakemake didn't print (it stops at the first missing one, which hides the scheduler log when a job died before writing its own) are read from disk
    - falls back to the previous message if nothing useful is found
- jobs that fail and are retried by snakemake (`retries:` in HPC profiles, or `--retries` in `--snakemake`) no longer stop the progress bar or make the workflow report a failure
  - a failure is only reported once snakemake gives up on the job. Previously the first error ended monitoring (and a workflow that went on to succeed was reported as failed)
  - only the error of the final attempt is passed to the error parser, instead of one for every attempt
- `harpy resume` no longer overwrites the harpy version of `workflow.yaml`
- error printing when using `--container` correctly displays full apptainer-prefixed shell call
- mitigated possibility of concurrent notebooks clashing when running on HPC

# Internal
- simplified summaries logic
- harpy-utils: replace Golang `regexp` with `coregex` for speed/efficiency
- progressbar logic moved out of `HarpyPrint` class in `printing.py` and into `progress.py`
  - removed the redundant `rich.Live` wrap to the progress bars-- progress bars should be a little snappier
- swapped pandas for polars (speed!)
- hidden command `harpy-utils process-notebook` replaced with `harpy-utils run-notebook`, which combines a nuanced
python-API call to `papermill` with the post-processing that was previously covered by `process-notebook`. 
  - functionally, this means the command line interface is cleaner in workflows, and the kernel engine can be declared to be IPC instead of TCP, which works better for concurrency.


# Added but not exposed yet
- `harpy snp deepvariant` - call SNPs in high-depth samples using Google's DeepVariant. AI AI AI!!!
  - this option must use the docker image of the software because of the way DeepVariant is packaged, so `--container` isn't exposed
- this is just internal scaffolding for now. Can be made public with sufficient interest
- `harpy sv genotype` - adds support for genotyping known sv breakpoints. 
  - this will be made public after a outstanding Issues/PRs are resolved