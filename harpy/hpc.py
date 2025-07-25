"""Harpy module to create HPC profiles for running snakemake on a cluster"""

import os
import subprocess
import sys
import rich_click as click
from rich import print as rprint
from rich.markdown import Markdown
from .common.printing import print_notice

def package_exists(pkg):
    """helper function to search for a package in the active conda environment"""
    search_pkg = subprocess.run(f"conda list snakemake-executor-plugin-{pkg}".split(), capture_output = True, text = True)
    exists = False
    for i in search_pkg.stdout.split("\n"):
        if i and not i.startswith("#"):
            exists = True
    if not exists:
        print_notice(Markdown(f"""
Using this scheduler requires installing a Snakemake plugin which wasn't detected in this environment.
It can be installed with:
    
```bash
conda install bioconda::snakemake-executor-plugin-{pkg}
# or #
pixi add snakemake-executor-plugin-{pkg}
```
    """))

@click.command()
def hpc_generic():
    """
    Create a template config for a generic scheduler
    
    This command creates a configuration for a generic HPC scheduler. Writes to `stdout`.
    You will also need to install `snakemake-executor-plugin-generic` for the HPC job submission to work.
    """
    sys.stdout.write("__use_yte__: true\n")
    sys.stdout.write("executor: cluster-generic\n")
    sys.stdout.write("default-resources:\n")
    sys.stdout.write("  mem_mb: attempt * 2000\n")
    sys.stdout.write("jobs: 50\n")
    sys.stdout.write("latency-wait: 60\n")
    sys.stdout.write("retries: 1\n")
    sys.stdout.write("# command for submitting jobs\n")
    sys.stdout.write("cluster-generic-submit-cmd: VALUE\n")
    sys.stdout.write("# [optional] command for retrieving job status\n")
    sys.stdout.write("cluster-generic-status-cmd: VALUE\n")
    sys.stdout.write("# [optional] command for cancelling jobs-- expected to take one or more jobids as arguments\n")
    sys.stdout.write("cluster-generic-cancel-cmd: VALUE\n")
    sys.stdout.write("# [optional] number of jobids to pass to cancel_cmd and if more are given, cancel_cmd will be called multiple times\n")
    sys.stdout.write("cluster-generic-cancel-nargs: 20\n")
    sys.stdout.write("# [optional] command for sidecar process\n")
    sys.stdout.write("cluster-generic-sidecar-cmd: VALUE\n")
    sys.stdout.write("\n# This section is for advanced copying into a scratch directory #\n")
    sys.stdout.write("## requires snakemake-storage-plugin-fs, which can be installed via conda\n")
    sys.stdout.write("#default-storage-provider: fs\n")
    sys.stdout.write("#local-storage-prefix: /home2/$USER\n")
    sys.stdout.write("#shared-fs-usage:\n")
    sys.stdout.write("#  - persistence\n")
    sys.stdout.write("#  - software-deployment\n")
    sys.stdout.write("#  - sources\n")
    sys.stdout.write("#  - source-cache\n")
    package_exists("cluster-generic")

@click.command()
def hpc_lsf():
    """
    Create a template config for LSF
    
    This command creates a configuration for the LSF HPC scheduler. Writes to `stdout`.
    You will also need to install `snakemake-executor-plugin-lsf` for the HPC job submission to work.
    """
    sys.stdout.write("__use_yte__: true\n")
    sys.stdout.write("executor: lsf\n")
    sys.stdout.write("default-resources:\n")
    sys.stdout.write("  lsf_queue: \n")
    sys.stdout.write("  walltime: 60 # minutes per job\n")
    sys.stdout.write("  mem_mb: attempt * 2000\n")
    sys.stdout.write("  # other args to pass to bsub\n")
    sys.stdout.write("  lsf_extra: VALUE\n")
    sys.stdout.write("jobs: 50\n")
    sys.stdout.write("latency-wait: 60\n")
    sys.stdout.write("retries: 1\n")
    sys.stdout.write("\n# This section is for advanced copying into a scratch directory #\n")
    sys.stdout.write("## requires snakemake-storage-plugin-fs, which can be installed via conda\n")
    sys.stdout.write("#default-storage-provider: fs\n")
    sys.stdout.write("#local-storage-prefix: /home2/$USER\n")
    sys.stdout.write("#shared-fs-usage:\n")
    sys.stdout.write("#  - persistence\n")
    sys.stdout.write("#  - software-deployment\n")
    sys.stdout.write("#  - sources\n")
    sys.stdout.write("#  - source-cache\n")
    package_exists("lsf")

@click.command()
def hpc_slurm():
    """
    Create a template config for SLURM
    
    This command creates a configuration for the SLURM HPC scheduler. Writes to `stdout`.
    You will also need to install `snakemake-executor-plugin-slurm` for the HPC job submission to work.
    """
    sys.stdout.write("__use_yte__: true\n")
    sys.stdout.write("executor: slurm\n")
    sys.stdout.write("default-resources:\n")
    sys.stdout.write("  slurm_account: $USER\n")
    sys.stdout.write("  slurm_partition: regular\n")
    sys.stdout.write("  mem_mb: attempt * 2000\n")
    sys.stdout.write("  runtime: 10\n")
    sys.stdout.write("jobs: 50\n")
    sys.stdout.write("latency-wait: 60\n")
    sys.stdout.write("retries: 1\n")
    sys.stdout.write("\n# This section is for advanced copying into a scratch directory #\n")
    sys.stdout.write("## requires snakemake-storage-plugin-fs, which can be installed via conda\n")
    sys.stdout.write("#default-storage-provider: fs\n")
    sys.stdout.write("#local-storage-prefix: /home2/$USER\n")
    sys.stdout.write("#shared-fs-usage:\n")
    sys.stdout.write("#  - persistence\n")
    sys.stdout.write("#  - software-deployment\n")
    sys.stdout.write("#  - sources\n")
    sys.stdout.write("#  - source-cache\n")
    package_exists("slurm")

@click.command()
def hpc_googlebatch():
    """
    Create a template config for Google Batch
    
    This command creates a configuration (`hpc/googlebatch.yaml`) for the Google Batch scheduler.
    You will also need to install `snakemake-executor-plugin-googlebatch` for the HPC job submission to work.
    """
    sys.stdout.write("__use_yte__: true\n")
    sys.stdout.write("executor: googlebatch\n")
    sys.stdout.write("jobs: 50\n")
    sys.stdout.write("latency-wait: 45\n")
    sys.stdout.write("retries: 1\n")
    sys.stdout.write("default-resources:\n")
    sys.stdout.write("## YOU MAY NOT NEED ALL OF THESE! ##\n")
    sys.stdout.write("# The name of the Google Project\n")
    sys.stdout.write("  googlebatch_project: Harpy\n")
    sys.stdout.write("# The name of the Google Project region (e.g., 'us-central1')\n")
    sys.stdout.write("  googlebatch_region: 'us-central1'\n")
    sys.stdout.write("# Retry count\n")
    sys.stdout.write("  googlebatch_retry_count: 1\n")
    sys.stdout.write("# Maximum run duration, string (e.g., '3600s')\n")
    sys.stdout.write("  googlebatch_max_run_duration: '3600s'\n")
    sys.stdout.write("# Memory in MiB\n")
    sys.stdout.write("  googlebatch_memory: attempt * 2000\n")
    sys.stdout.write("# The default number of work tasks (these are NOT MPI ranks)\n")
    sys.stdout.write("  googlebatch_work_tasks: 50\n")
    sys.stdout.write("# The default number of work tasks per node (NOT MPI ranks)\n")
    sys.stdout.write("  googlebatch_work_tasks_per_node: 10\n")
    sys.stdout.write("# Milliseconds per cpu-second\n")
    sys.stdout.write("  googlebatch_cpu_milli: 1000\n")
    sys.stdout.write("# A custom container for use with Google Batch COS\n")
    sys.stdout.write("  googlebatch_container: VALUE\n")
    sys.stdout.write("# A docker registry password for COS if credentials are required\n")
    sys.stdout.write("  googlebatch_docker_password: VALUE\n")
    sys.stdout.write("# A docker registry username for COS if credentials are required\n")
    sys.stdout.write("  googlebatch_docker_username: VALUE\n")
    sys.stdout.write("# Google Cloud machine type or VM (mpitune on c2 and c2d family)\n")
    sys.stdout.write("  googlebatch_machine_type: 'c2-standard-4'\n")
    sys.stdout.write("# Comma separated key value pairs to label job (e.g., model=a3,stage=test)\n")
    sys.stdout.write("  googlebatch_labels: VALUE\n")
    sys.stdout.write("# Google Cloud image family (defaults to hpc-centos-7)\n")
    sys.stdout.write("  googlebatch_image_family: 'hpc-centos-7'\n")
    sys.stdout.write("# Selected image project\n")
    sys.stdout.write("  googlebatch_image_project: 'cloud-hpc-image-public'\n")
    sys.stdout.write("# Boot disk size (GB)\n")
    sys.stdout.write("  googlebatch_boot_disk_gb: VALUE\n")
    sys.stdout.write("# The URL of an existing network resource\n")
    sys.stdout.write("  googlebatch_network: VALUE\n")
    sys.stdout.write("# The URL of an existing subnetwork resource\n")
    sys.stdout.write("  googlebatch_subnetwork: VALUE\n")
    sys.stdout.write("# Boot disk type. (e.g., gcloud compute disk-types list)\n")
    sys.stdout.write("  googlebatch_boot_disk_type: VALUE\n")
    sys.stdout.write("# Boot disk image (e.g., batch-debian, bath-centos)\n")
    sys.stdout.write("  googlebatch_boot_disk_image: VALUE\n")
    sys.stdout.write("# Mount path for Google bucket (if defined)\n")
    sys.stdout.write("  googlebatch_mount_path: '/mnt/share'\n")
    sys.stdout.write("# One or more snippets to add to the Google Batch task setup\n")
    sys.stdout.write("  googlebatch_snippets: VALUE\n")
    package_exists("googlebatch")
