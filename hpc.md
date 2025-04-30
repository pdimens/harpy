---
label: Running on an HPC
icon: server
order: 5
---
# :icon-server: Running on an HPC Cluster

You'll be working with genomic data, so it's very likely you'll try to run Harpy on a high-performance computing (HPC) cluster at some point.
Doing so will usually require you to log into a head node and submit jobs that some kind of scheduler (e.g. SLURM, HTCondor)
will manage by sending to the worker nodes. Naturally, it would seem like running Harpy on an HPC should be done in the typical
way of writing a job script where you call Harpy and submitting that for execution, and there's a decent chance that might work
without any fuss.

However tempting (or functional) that may be, that is not the idomatic way of running a Harpy workflow on an HPC, and there may
be situations where a simple job script submission won't work, like if the worker nodes do not have internet access, or have specific 
network-mounted filesystem configurations (see [this discussion](https://github.com/pdimens/harpy/discussions/222#discussion-8113022)). All 
major Harpy workflows rely on Snakemake and, conveniently/thankfully, the Snakemake developers put a lot of effort into addressing workflow
execution on HPC clusters.

## Snakemake HPC execution
Snakemake has [executor plugins](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) for various cluster job 
managers (like SLURM) and storage plugins (like `fs` or `s3`) to automate workflow execution and file transfer in HPC contexts. In order
to activate the "HPC mode" (not a real phrase) of Snakemake, you need to provide a yaml file with configuration details for your HPC job
manager and [possibly] the file storage system. Behind the scenes, Snakemake will use the configuration information and 
**automatically create and submit** individual job scripts for **every single job** of the workflow on your behalf. It's kind of amazing.
As an example, a configuration could look something like this:

```yaml
executor: slurm
default-resources:
    slurm_account: "accountname"
    mem_mb_per_cpu: 1800
    runtime: "90m"

latency-wait: 5
default-storage-provider: fs
shared-fs-usage:
  - persistence
  - software-deployment
  - sources
  - source-cache
remote-job-local-storage-prefix: "/home2/accountname/SCRATCH"
local-storage-prefix: "/home/accountname/DATA"
```

## HPC features in Harpy
Harpy provides this Snakemake-driven HPC support with the `--hpc` option available to
most workflows. This option requires a path to the configuration file with the HPC
configuration yaml. In practice, that would look something like:

```bash
harpy qc -a auto --hpc hpc/slurm.yaml data/porcupine
```

In addition to the config file, you will need to install the executor plugins you intend to use. This is done with
e.g. `conda install bioconda::snakemake-executor-plugin-slurm` and ` conda install bioconda::snakemake-storage-plugin-fs ` or their
Pixi equivalents with e.g. `pixi add snakemake-executor-plugin-slurm`.

### Configuration templates
This configuration stuff is a lot of congitive burden in addition to just trying to process your data, so you can use
[!badge corners="pill" text="harpy template hpc-"](/Workflows/other.md/#hpc-)
to create skeleton configurations for various supported cluster managers and fill in the information you need:
```bash
harpy template hpc-lsf > lsf.config.yml
```

Depending on your system, it may be necessary to read the [documentation](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) 
for a particular executor plugin and understand what configuration options their API exposes. The configurations can start to become
very technical, so we recommend starting with a simple configuration and getting more complex if issues arise. The Snakemake
executor plugins are admittedly not consistent in their documentation quality and it's sometimes a rapidly
changing landscape (for example, the `HTCondor` plugin was recently and quietly deperacated). If an executor plugin exists
that you would like `harpy template` support for, please [open an Issue](https://github.com/pdimens/harpy/issues/new?template=feature_request.yml) and we'll get it added!