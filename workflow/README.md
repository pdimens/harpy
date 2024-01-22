# Snakemake Workflow

This directory follows the standard Snakemake directory structure 
in a way that makes for how Harpy runs as a command-line software.

### `envs`
Contains YAML files that are used during development (CI, conda build, containers)

### `report`
Contains the RMarkdown files used to create nice FlexDashboards by workflows

### `rules`
The snakefiles defining the different workflows

### `scripts`
The various custom scripts used by different rules