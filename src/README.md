# The Harpy package

This directory contains the Harpy python package. This contains the command-line
interface users interact with when they call `harpy qc ...`. It's where
file inputs are checked, validated, and which workflows are executed by snakemake.

This part of Harpy only has a dependency on `snakemake`, `rich-click` (and by
extension, `rich` and `click`), samtools, and bcftools.

## Contents
- `bin`: scripts that are not expected to be run inside a container and copied to the conda PATH
- `reports`: RMarkdown files that are used to create the nice flexdashboard reports
- `scripts`: scripts with dependencies that require a conda/container to have its dependencies available. these get copied into the `workdir/scripts` folder of a workflow
- `snakefiles`: the snakefiles that govern the workflows for most harpy modules
- everything else: the components of the harpy command line interface 