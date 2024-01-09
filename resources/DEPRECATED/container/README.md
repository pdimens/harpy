# Docker container configurations

This directory contains the environment information required for creating
Docker images preconfigured to work for the various Harpy workflows. Each
directory (`qc`, `align`, etc.) corresponds to that Harpy module (e.g. `qc/`
corresponds to the `harpy qc` module). Each of these folders contain three
files:
1. a `Dockerfile` that is used to create the container that lives on DockerHub
and will be pulled by Snakemake at runtime.
2. a `.yaml` file the specifies the conda environment (channels, dependencies)
required for the workflow that is used to generate the lock file `mamba` uses
inside the Docker image to quickly install the programs without having to
spend time solving the dependency graph each time.
3. a `.lock` file that is generated from the `.yaml` file by `conda-lock`

The combination of these things accomplishes:
1. A user will only need to install the most basic necessities when installing 
Harpy, which are essentially Harpy, Snakemake, and Rich-Click.
2. The remaining dependencies will be offloaded into containers that are more
HPC-admin friendly and will save users a lot of installation time and 
overall disk space. 
3. Dependencies are lumped by workflow, meaning it will be easier to update 
individual workflows instead of relying on a global set of packages that work 
for all workflows. And while yes you can do this with Snakemake + conda environments, it will create a lot more redundancy and disk usage.
4. The yaml-lockfile combination shaves off a lot of time in Docker image creation, especially during development.