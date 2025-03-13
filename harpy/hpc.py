"""Harpy module to create HPC profiles for running snakemake on a cluster"""

import os
import subprocess
import sys
import rich_click as click
from rich import print as rprint
from rich.markdown import Markdown
from ._printing import print_notice

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def hpc():
    """
    Profile templates for cluster job submissions.
    
    This command group provides subcommands to generate scheduler-specific YAML
    configuration files (e.g., "generic.yaml", "lsf.yaml", "htcondor.yaml", "slurm.yaml",
    or "googlebatch.yaml") in the "hpc/" directory. These profiles enable Snakemake to manage
    job submissions on HPC clusters. Ensure that the corresponding Snakemake executor plugin
    is installed for the intended scheduler.
    """

docstring = {
    "harpy hpc": [
        {
            "name": "Submission Systems",
            "commands": ["generic", "googlebatch", "htcondor", "lsf", "slurm"],
            "panel_styles": {"border_style": "blue"}
        },
    ]
}

def package_exists(pkg):
    """Verify the installation of a required Snakemake executor plugin.
    
    This function checks if a package named "snakemake-executor-plugin-{pkg}" is
    installed in the active conda environment. If not, it prints a markdown-formatted
    notice with instructions to install the missing plugin using either conda or pixi.
    
    Args:
        pkg: Suffix of the plugin package name; the complete name is formed as
             "snakemake-executor-plugin-{pkg}".
    """
    search_pkg = subprocess.run(f"conda list snakemake-executor-plugin-{pkg}".split(), capture_output = True, text = True)
    exists = False
    for i in search_pkg.stdout.split("\n"):
        if i and not i.startswith("#"):
            exists = True
    if not exists:
        print_notice(Markdown(f"""
Using this scheduler requires installing an additional Snakemake plugin, which wasn't detected in this environment.
Install the `{pkg}` plugin with:
    
```bash
conda install bioconda::snakemake-executor-plugin-{pkg}
# or #
pixi add snakemake-executor-plugin-{pkg}
```
    """))

@click.command()
def generic():
    """
    Generates a YAML configuration file for a generic scheduler.
    
    This function creates or overwrites the "hpc/generic.yaml" file with default
    configuration settings used by the generic HPC scheduler for Snakemake. It ensures
    the necessary "hpc" directory exists, verifies that the "cluster-generic" package is
    installed, and writes key settings such as memory allocation, job limits, and commands
    for job submission, status querying, cancellation, and sidecar processing.
    """
    outfile = "hpc/generic.yaml"
    os.makedirs("hpc", exist_ok=True)
    if os.path.exists(outfile):
        rprint(f"HPC profile [blue]{outfile}[/blue] already exists, overwriting\n", file = sys.stderr)
    package_exists("cluster-generic")
    with open(outfile, "w", encoding = "utf-8") as yml:
        yml.write("__use_yte__: true\n")
        yml.write("executor: cluster-generic\n")
        yml.write("default-resources:\n")
        yml.write("  mem_mb: attempt * 2000\n")
        yml.write("jobs: 50\n")
        yml.write("latency-wait: 60\n")
        yml.write("retries: 1\n")
        yml.write("# command for submitting jobs\n")
        yml.write("cluster-generic-submit-cmd: VALUE\n")
        yml.write("# [optional] command for retrieving job status\n")
        yml.write("cluster-generic-status-cmd: VALUE\n")
        yml.write("# [optional] command for cancelling jobs-- expected to take one or more jobids as arguments\n")
        yml.write("cluster-generic-cancel-cmd: VALUE\n")
        yml.write("# [optional] number of jobids to pass to cancel_cmd and if more are given, cancel_cmd will be called multiple times\n")
        yml.write("cluster-generic-cancel-nargs: 20\n")
        yml.write("# [optional] command for sidecar process\n")
        yml.write("cluster-generic-sidecar-cmd: VALUE\n")
        yml.write("\n# This section is for advanced copying into a scratch directory #\n")
        yml.write("## requires snakemake-storage-plugin-fs, which can be installed via conda\n")
        yml.write("#default-storage-provider: fs\n")
        yml.write("#local-storage-prefix: /home2/$USER\n")
        yml.write("#shared-fs-usage:\n")
        yml.write("#  - persistence\n")
        yml.write("#  - software-deployment\n")
        yml.write("#  - sources\n")
        yml.write("#  - source-cache\n")

@click.command()
def lsf():
    """
    Generate a YAML configuration file for the LSF scheduler.
    
    This function creates the 'hpc' directory (if it does not exist) and writes the LSF
    configuration to 'hpc/lsf.yaml'. It verifies that the required LSF executor package is
    installed and overwrites any existing configuration file. The output YAML includes default
    resource settings and job submission parameters, with placeholder values for further customization.
    """
    os.makedirs("hpc", exist_ok=True)
    outfile = "hpc/lsf.yaml"
    if os.path.exists(outfile):
        rprint(f"HPC profile [blue]{outfile}[/blue] already exists, overwriting\n", file = sys.stderr)
    package_exists("lsf")
    with open(outfile, "w", encoding = "utf-8") as yml:
        yml.write("__use_yte__: true\n")
        yml.write("executor: lsf\n")
        yml.write("default-resources:\n")
        yml.write("  lsf_queue: \n")
        yml.write("  walltime: 60 # minutes per job\n")
        yml.write("  mem_mb: attempt * 2000\n")
        yml.write("  # other args to pass to bsub\n")
        yml.write("  lsf_extra: VALUE\n")
        yml.write("jobs: 50\n")
        yml.write("latency-wait: 60\n")
        yml.write("retries: 1\n")
        yml.write("\n# This section is for advanced copying into a scratch directory #\n")
        yml.write("## requires snakemake-storage-plugin-fs, which can be installed via conda\n")
        yml.write("#default-storage-provider: fs\n")
        yml.write("#local-storage-prefix: /home2/$USER\n")
        yml.write("#shared-fs-usage:\n")
        yml.write("#  - persistence\n")
        yml.write("#  - software-deployment\n")
        yml.write("#  - sources\n")
        yml.write("#  - source-cache\n")

@click.command()
def htcondor():
    """
    Generate a YAML configuration file for the HTCondor scheduler.
    
    This function ensures the "hpc" directory exists and creates or overwrites the 
    "hpc/htcondor.yaml" file with default configuration settings for running Snakemake 
    jobs on an HTCondor cluster. It warns if an existing file is being overwritten and 
    checks that the HTCondor executor plugin is installed. The configuration includes 
    default resource parameters and a commented section for advanced scratch directory setup.
    """
    os.makedirs("hpc", exist_ok=True)
    outfile = "hpc/htcondor.yaml"
    if os.path.exists(outfile):
        rprint(f"HPC profile [blue]{outfile}[/blue] already exists, overwriting\n", file = sys.stderr)
    package_exists("htcondor")
    with open(outfile, "w", encoding = "utf-8") as yml:
        yml.write("__use_yte__: true\n")
        yml.write("executor: htcondor\n")
        yml.write("default-resources:\n")
        yml.write("  getenv: True\n")
        yml.write("  rank: \n")
        yml.write("  max-retries: 1\n")
        yml.write("  request_memory: attempt * 2000\n")
        yml.write("jobs: 50\n")
        yml.write("latency-wait: 60\n")
        yml.write("retries: 1\n")
        yml.write("\n# This section is for advanced copying into a scratch directory #\n")
        yml.write("## requires snakemake-storage-plugin-fs, which can be installed via conda\n")
        yml.write("#default-storage-provider: fs\n")
        yml.write("#local-storage-prefix: /home2/$USER\n")
        yml.write("#shared-fs-usage:\n")
        yml.write("#  - persistence\n")
        yml.write("#  - software-deployment\n")
        yml.write("#  - sources\n")
        yml.write("#  - source-cache\n")

@click.command()
def slurm():
    """
    Generate a SLURM configuration file for HPC execution.
    
    This function creates or overwrites the "hpc/slurm.yaml" file with default
    configuration settings for running jobs on a SLURM scheduler. It ensures the
    required SLURM plugin is installed and writes YAML-formatted parameters,
    including default resource allocations, job count, latency wait, and retry
    settings. Commented sections provide options for configuring advanced scratch
    directory usage.
    """
    os.makedirs("hpc", exist_ok=True)
    outfile = "hpc/slurm.yaml"
    if os.path.exists(outfile):
        rprint(f"HPC profile [blue]{outfile}[/blue] already exists, overwriting\n", file = sys.stderr)
    package_exists("slurm")
    with open(outfile, "w", encoding = "utf-8") as yml:
        yml.write("__use_yte__: true\n")
        yml.write("executor: slurm\n")
        yml.write("default-resources:\n")
        yml.write("  slurm_account: $USER\n")
        yml.write("  slurm_partition: regular\n")
        yml.write("  mem_mb: attempt * 2000\n")
        yml.write("  runtime: 10\n")
        yml.write("jobs: 50\n")
        yml.write("latency-wait: 60\n")
        yml.write("retries: 1\n")
        yml.write("\n# This section is for advanced copying into a scratch directory #\n")
        yml.write("## requires snakemake-storage-plugin-fs, which can be installed via conda\n")
        yml.write("#default-storage-provider: fs\n")
        yml.write("#local-storage-prefix: /home2/$USER\n")
        yml.write("#shared-fs-usage:\n")
        yml.write("#  - persistence\n")
        yml.write("#  - software-deployment\n")
        yml.write("#  - sources\n")
        yml.write("#  - source-cache\n")

@click.command()
def googlebatch():
    """Generate a Google Batch HPC profile configuration file.
    
    Creates or overwrites the 'hpc/googlebatch.yaml' file with default settings for
    Google Batch job submission. Verifies that the required 'googlebatch' plugin is
    installed before writing the configuration.
    """
    os.makedirs("hpc", exist_ok=True)
    outfile = "hpc/googlebatch.yaml"
    if os.path.exists(outfile):
        rprint(f"HPC profile [blue]{outfile}[/blue] already exists, overwriting\n", file = sys.stderr)
    package_exists("googlebatch")
    with open(outfile, "w", encoding = "utf-8") as yml:
        yml.write("__use_yte__: true\n")
        yml.write("executor: googlebatch\n")
        yml.write("jobs: 50\n")
        yml.write("latency-wait: 45\n")
        yml.write("retries: 1\n")
        yml.write("default-resources:\n")
        yml.write("## YOU MAY NOT NEED ALL OF THESE! ##\n")
        yml.write("# The name of the Google Project\n")
        yml.write("  googlebatch_project: Harpy\n")
        yml.write("# The name of the Google Project region (e.g., 'us-central1')\n")
        yml.write("  googlebatch_region: 'us-central1'\n")
        yml.write("# Retry count\n")
        yml.write("  googlebatch_retry_count: 1\n")
        yml.write("# Maximum run duration, string (e.g., '3600s')\n")
        yml.write("  googlebatch_max_run_duration: '3600s'\n")
        yml.write("# Memory in MiB\n")
        yml.write("  googlebatch_memory: attempt * 2000\n")
        yml.write("# The default number of work tasks (these are NOT MPI ranks)\n")
        yml.write("  googlebatch_work_tasks: 50\n")
        yml.write("# The default number of work tasks per node (NOT MPI ranks)\n")
        yml.write("  googlebatch_work_tasks_per_node: 10\n")
        yml.write("# Milliseconds per cpu-second\n")
        yml.write("  googlebatch_cpu_milli: 1000\n")
        yml.write("# A custom container for use with Google Batch COS\n")
        yml.write("  googlebatch_container: VALUE\n")
        yml.write("# A docker registry password for COS if credentials are required\n")
        yml.write("  googlebatch_docker_password: VALUE\n")
        yml.write("# A docker registry username for COS if credentials are required\n")
        yml.write("  googlebatch_docker_username: VALUE\n")
        yml.write("# Google Cloud machine type or VM (mpitune on c2 and c2d family)\n")
        yml.write("  googlebatch_machine_type: 'c2-standard-4'\n")
        yml.write("# Comma separated key value pairs to label job (e.g., model=a3,stage=test)\n")
        yml.write("  googlebatch_labels: VALUE\n")
        yml.write("# Google Cloud image family (defaults to hpc-centos-7)\n")
        yml.write("  googlebatch_image_family: 'hpc-centos-7'\n")
        yml.write("# Selected image project\n")
        yml.write("  googlebatch_image_project: 'cloud-hpc-image-public'\n")
        yml.write("# Boot disk size (GB)\n")
        yml.write("  googlebatch_boot_disk_gb: VALUE\n")
        yml.write("# The URL of an existing network resource\n")
        yml.write("  googlebatch_network: VALUE\n")
        yml.write("# The URL of an existing subnetwork resource\n")
        yml.write("  googlebatch_subnetwork: VALUE\n")
        yml.write("# Boot disk type. (e.g., gcloud compute disk-types list)\n")
        yml.write("  googlebatch_boot_disk_type: VALUE\n")
        yml.write("# Boot disk image (e.g., batch-debian, bath-centos)\n")
        yml.write("  googlebatch_boot_disk_image: VALUE\n")
        yml.write("# Mount path for Google bucket (if defined)\n")
        yml.write("  googlebatch_mount_path: '/mnt/share'\n")
        yml.write("# One or more snippets to add to the Google Batch task setup\n")
        yml.write("  googlebatch_snippets: VALUE\n")

hpc.add_command(slurm)
hpc.add_command(htcondor)
hpc.add_command(lsf)
hpc.add_command(generic)
hpc.add_command(googlebatch)
