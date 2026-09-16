# New
- `harpy view envs` prints versions too, where applicable
- CLI validations for `--hpc`/`-H` to check if the plugins necessary for the configuration are installed
- `harpy preprocess meier2021` now uses `dmox` v0.3.1
  - introduces the `--use-stitch-base` logic via `--stitch`
  - introduces the `--allow-complementary-stitch` logic via `--stitch-comp`
- `harpy snp deepvariant` - call SNPs in high-depth samples using Google's DeepVariant. AI AI AI!!!
  - this option must use the docker image of the software because of the way DeepVariant is packaged, so `--container` isn't exposed

# Fixes
- `harpy view envs`: simpler logic and print diagnostic text if empty
- more robust snakemake error printing (again)
  - this time it's a parse-and-gather approach that uses an internal class
  - it should print scheduler messages on error now (maybe? you never know with snakemake)
  - it should print resources on error now (again, maybe?)
  - strengthened outputting snakemake missing and syntax errors
- `harpy resume` no longer overwrites the harpy version of `workflow.yaml`
- error printing when using `--container` correctly displays full apptainer invocation

# Breaking
- the `dmox` update makes `workflow.yaml` files from previous harpy versions invalid
  - not technically invalid, but the key absence/mismatch will default to all optional features turned off, which may be unintended

# Internal
- simplifies summaries logic
- prep arachne workflow
- harpy-utils: replace Golang `regexp` with `coregex` for speed/efficiency
- progressbar logic moved out of `HarpyPrint` class in `printing.py` and into `progress.py`
  - removed the redundant `rich.Live` wrap to the progress bars-- progress bars should be a little snappier