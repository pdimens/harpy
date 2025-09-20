"""Harpy module to create HPC profiles for running snakemake on a cluster"""

import sys
import rich_click as click
from harpy.common.system_ops import package_absent

@click.command(panel = "HPC Configurations")
def hpc_generic():
    """
    Create a template config for a generic scheduler
    
    This command creates a configuration for a generic HPC scheduler. Writes to `stdout`.
    You will also need to install `snakemake-executor-plugin-cluster-generic` for the HPC job submission to work.
    """
    out_text = []
    out_text.append("__use_yte__: true")
    out_text.append("executor: cluster-generic")
    out_text.append("default-resources:")
    out_text.append("  mem_mb: attempt * 3200")
    out_text.append('  time: "12:00:00"')
    out_text.append('  mem_per_cpu: attempt * 3200')
    out_text.append('  tmpdir: "/tmp"')
    out_text.append("jobs: 100")
    out_text.append("max-jobs-per-second: 10")
    out_text.append("latency-wait: 60")
    out_text.append("retries: 1")
    out_text.append("\n# command for submitting jobs")
    out_text.append("cluster-generic-submit-cmd: VALUE")
    out_text.append("\n# [optional] command for retrieving job status")
    out_text.append("cluster-generic-status-cmd: VALUE")
    out_text.append("\n# [optional] command for cancelling jobs-- expected to take one or more jobids as arguments")
    out_text.append("cluster-generic-cancel-cmd: VALUE")
    out_text.append("\n# [optional] number of jobids to pass to cancel_cmd and if more are given, cancel_cmd will be called multiple times")
    out_text.append("cluster-generic-cancel-nargs: 20")
    out_text.append("\n# [optional] command for sidecar process")
    out_text.append("cluster-generic-sidecar-cmd: VALUE")
    out_text.append("\n# This section is for advanced copying into a scratch directory #")
    out_text.append("## requires snakemake-storage-plugin-fs, which can be installed via conda")
    out_text.append("#default-storage-provider: fs")
    out_text.append("#local-storage-prefix: /home2/$USER")
    out_text.append("#shared-fs-usage:")
    out_text.append("#  - persistence")
    out_text.append("#  - software-deployment")
    out_text.append("#  - sources")
    out_text.append("#  - source-cache")
    sys.stdout.write("\n".join(out_text) + "\n")

    package_absent("snakemake-executor-plugin-cluster-generic")

#def hpc_generic2():
#    pass
    #executor: cluster-generic
    #cluster-generic-submit-cmd:
    #  mkdir -p results/slurm_logs/{rule} &&
    #  sbatch
    #    --partition=regular
    #    --account=username
    #    --cpus-per-task={threads}
    #    --mem-per-cpu={resources.mem_per_cpu}
    #    --time={resources.time}
    #    --job-name=harpy-{rule}-{wildcards}
    #    --output=results/slurm_logs/{rule}/{rule}-%j-{wildcards}.out
    #    --error=results/slurm_logs/{rule}/{rule}-%j-{wildcards}.err
    #    --parsable
    #cluster-generic-status-cmd: status-sacct-robust.sh
    #cluster-generic-cancel-cmd: scancel
    #cluster-generic-cancel-nargs: 400
    #default-resources:
    #  - time="12:00:00"
    #  - mem_per_cpu=3200
    #  - tmpdir="/tmp"
    #restart-times: 2
    #max-jobs-per-second: 10
    #max-status-checks-per-second: 2
    #local-cores: 1
    #latency-wait: 60
    #cores: 800
    #jobs: 500
    #keep-going: True
    #rerun-incomplete: True

@click.command(panel = "HPC Configurations")
def hpc_lsf():
    """
    Create a template config for LSF
    
    This command creates a configuration for the LSF HPC scheduler. Writes to `stdout`.
    You will also need to install `snakemake-executor-plugin-lsf` for the HPC job submission to work.
    """
    out_text = []
    out_text.append("__use_yte__: true")
    out_text.append("executor: lsf")
    out_text.append("default-resources:")
    out_text.append("  lsf_queue: ")
    out_text.append("  walltime: 60 # minutes per job")
    out_text.append("  mem_mb: attempt * 2000")
    out_text.append("  # other args to pass to bsub")
    out_text.append("  lsf_extra: VALUE")
    out_text.append("jobs: 50")
    out_text.append("latency-wait: 60")
    out_text.append("retries: 1")
    out_text.append("\n# This section is for advanced copying into a scratch directory #")
    out_text.append("## requires snakemake-storage-plugin-fs, which can be installed via conda")
    out_text.append("#default-storage-provider: fs")
    out_text.append("#local-storage-prefix: /home2/$USER")
    out_text.append("#shared-fs-usage:")
    out_text.append("#  - persistence")
    out_text.append("#  - software-deployment")
    out_text.append("#  - sources")
    out_text.append("#  - source-cache")
    sys.stdout.write("\n".join(out_text) + "\n")

    package_absent("snakemake-executor-plugin-lsf")

@click.command(panel = "HPC Configurations")
def hpc_slurm():
    """
    Create a template config for SLURM
    
    This command creates a configuration for the SLURM HPC scheduler. Writes to `stdout`.
    You will also need to install `snakemake-executor-plugin-slurm` for the HPC job submission to work.
    """
    out_text = []
    out_text.append("__use_yte__: true")
    out_text.append("executor: slurm")
    out_text.append("default-resources:")
    out_text.append("  slurm_account: $USER")
    out_text.append("  slurm_partition: regular")
    out_text.append("  mem_mb: attempt * 2000")
    out_text.append("  runtime: 10")
    out_text.append("jobs: 50")
    out_text.append("latency-wait: 60")
    out_text.append("retries: 1")
    out_text.append("\n# This section is for advanced copying into a scratch directory #")
    out_text.append("## requires snakemake-storage-plugin-fs, which can be installed via conda")
    out_text.append("#default-storage-provider: fs")
    out_text.append("#local-storage-prefix: /home2/$USER")
    out_text.append("#shared-fs-usage:")
    out_text.append("#  - persistence")
    out_text.append("#  - software-deployment")
    out_text.append("#  - sources")
    out_text.append("#  - source-cache")
    sys.stdout.write("\n".join(out_text) + "\n")
    package_absent("snakemake-executor-plugin-slurm")

@click.command(panel = "HPC Configurations")
def hpc_googlebatch():
    """
    Create a template config for Google Batch
    
    This command creates a configuration for the Google Batch scheduler. Writes to `stdout`.
    You will also need to install `snakemake-executor-plugin-googlebatch` for the HPC job submission to work.
    """
    out_text = []
    out_text.append("__use_yte__: true")
    out_text.append("executor: googlebatch")
    out_text.append("jobs: 50")
    out_text.append("latency-wait: 45")
    out_text.append("retries: 1")
    out_text.append("default-resources:")
    out_text.append("## YOU MAY NOT NEED ALL OF THESE! ##")
    out_text.append("# The name of the Google Project")
    out_text.append("  googlebatch_project: Harpy")
    out_text.append("\n# The name of the Google Project region (e.g., 'us-central1')")
    out_text.append("  googlebatch_region: 'us-central1'")
    out_text.append("\n# Retry count")
    out_text.append("  googlebatch_retry_count: 1")
    out_text.append("\n# Maximum run duration, string (e.g., '3600s')")
    out_text.append("  googlebatch_max_run_duration: '3600s'")
    out_text.append("\n# Memory in MiB")
    out_text.append("  googlebatch_memory: attempt * 2000")
    out_text.append("\n# The default number of work tasks (these are NOT MPI ranks)")
    out_text.append("  googlebatch_work_tasks: 50")
    out_text.append("\n# The default number of work tasks per node (NOT MPI ranks)")
    out_text.append("  googlebatch_work_tasks_per_node: 10")
    out_text.append("\n# Milliseconds per cpu-second")
    out_text.append("  googlebatch_cpu_milli: 1000")
    out_text.append("\n# A custom container for use with Google Batch COS")
    out_text.append("  googlebatch_container: VALUE")
    out_text.append("\n# A docker registry password for COS if credentials are required")
    out_text.append("  googlebatch_docker_password: VALUE")
    out_text.append("\n# A docker registry username for COS if credentials are required")
    out_text.append("  googlebatch_docker_username: VALUE")
    out_text.append("\n# Google Cloud machine type or VM (mpitune on c2 and c2d family)")
    out_text.append("  googlebatch_machine_type: 'c2-standard-4'")
    out_text.append("\n# Comma separated key value pairs to label job (e.g., model=a3,stage=test)")
    out_text.append("  googlebatch_labels: VALUE")
    out_text.append("\n# Google Cloud image family (defaults to hpc-centos-7)")
    out_text.append("  googlebatch_image_family: 'hpc-centos-7'")
    out_text.append("\n# Selected image project")
    out_text.append("  googlebatch_image_project: 'cloud-hpc-image-public'")
    out_text.append("\n# Boot disk size (GB)")
    out_text.append("  googlebatch_boot_disk_gb: VALUE")
    out_text.append("\n# The URL of an existing network resource")
    out_text.append("  googlebatch_network: VALUE")
    out_text.append("\n# The URL of an existing subnetwork resource")
    out_text.append("  googlebatch_subnetwork: VALUE")
    out_text.append("\n# Boot disk type. (e.g., gcloud compute disk-types list)")
    out_text.append("  googlebatch_boot_disk_type: VALUE")
    out_text.append("\n# Boot disk image (e.g., batch-debian, bath-centos)")
    out_text.append("  googlebatch_boot_disk_image: VALUE")
    out_text.append("\n# Mount path for Google bucket (if defined)")
    out_text.append("  googlebatch_mount_path: '/mnt/share'")
    out_text.append("\n# One or more snippets to add to the Google Batch task setup")
    out_text.append("  googlebatch_snippets: VALUE")
    sys.stdout.write("\n".join(out_text) + "\n")

    package_absent("snakemake-executor-plugin-googlebatch")
