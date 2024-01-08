# The Harpy package

This directory contains the Harpy python package. This contains the command-line
interface users interact with when they call `harpy qc ...`. It's where
file inputs are checked, validated, and which workflows are executed by snakemake.

This "part" of Harpy only has a dependency on `snakemake` and `rich-click` (and by
extension, `rich` and `click`).