# New
- `harpy view envs`: print versions
- CLI validations for `--hpc`/`-H` to check if the plugins necessary for the configuration are installed

# Fixes
- `harpy view envs`: simpler logic and print diagnostic text if empty
- more robust snakemake error printing (again)
  - this time it's a parse-and-gather approach that uses an internal class
  - it should print scheduler messages on error now
  - it should print resources on error now
  - strengthened outputting snakemake missing and syntax errors
- `harpy resume` no longer overwrites the harpy version of `workflow.yaml`
- error printing when using `--container` correctly displays full apptainer invocation

# internal
- simplifies summaries logic
- prep arachne workflow
- replace Golang regexp with coregex for speed/efficiency