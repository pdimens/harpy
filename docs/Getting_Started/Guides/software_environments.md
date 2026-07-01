---
author: 
  - name: Pavel Dimens
    avatar: https://cals.cornell.edu/sites/default/files/styles/faculty/public/2024-09/afs-headshot-high-res-2cropped_0.jpg
date: 2024-06-21
category: guides
description: Deciding between using Conda or Containers
icon: container
label: Software Management
image: https://visualpharm.com/assets/917/Docker-595b40b65ba036ed117d3f62.svg
---

# :icon-container: Choosing a software runtime method
## TL;DR
- container is more likely to work on all systems, but much slower
  - `Apptainer` (what Snakemake uses) is only available on Linux systems
- conda is quicker and better for troubleshooting, but may have unexpected errors
  - but works on any operation system, assuming it's installed

There are two ways you can run Harpy, using a container with the necessary
software environments in it (`--container`), or with local conda environments
(the default). If software development and containerization 
isn't your jam, that's great, you're in the right place! Below is a quick
explanation of what/why and the tradeoffs between either approach so you
can decide for yourself which makes more sense to use.


## What Harpy Provides
An conda-based installation of Harpy provides only the minimal set of 
programs Harpy needs to begin a workflow. These include: `python`, `snakemake-minimal`, the `htslib` programs (e.g. `samtools`), along
with a few other smaller packages. Noticeably, there aren't sequence aligners, quality-assessment tools, phasers, etc. This is because 
some of the software dependencies themselves have clashing dependencies and cannot be installed alongside each other, but more 
importantly, it keeps the Harpy installation quite small and quick.

## How Harpy Provides the Other Stuff
Instead of a monolithic Harpy environment, which would be impossible with 
the current software dependencies, there are a handful of defined conda environment recipes that Harpy
workflows generate. Snakemake will make environments of those recipes, then jump in and out of those local conda 
environments as dictated by the software needs of any given job (given in 
the `conda:` directive  within a rule). Those local environments live inside 
`.environments/wildhashnumber`, with auto-generated names reflecting the 
hash of the environment (e.g. `.environments/21ceb8c2fe7dd21206ab90c2af8f847f_`).

**But**, those environments need to be created at runtime if they don't 
already exist in `.environments/`, so Harpy (technically Snakemake) will install 
them before running the jobs within a workflow. On some HPC systems, this
process can move glacially slow (it might be a RAID or NAS thing) and this
might make you think a Harpy workflow is hanging at the environment 
installation step before it even begins its first job. That isn't ideal.
Additionally, sysadmins aren't particularly fond of how many files are 
created with conda-based installations, which leads us to containerization.
**Note**: Pixi tends to be significantly faster for this kind of thing and when
Snakemake officially integrates Pixi into workflow management, Harpy will make the
shift shortly after.

## Harpy and Containers
==- Containers, a Primer
If you aren't sure exactly what containers are, great, we aren't either! But
here's what we do know: it's a tiny mountable file containing an entire
operating system and whatever other bits you might need. Creating containers
is done with a recipe that takes a base "image" (an established existing 
container) and adds "layers" of modifications to that base image. Imagine a 
simple recipe where you declare a base image of a minimal Ubuntu 22 system 
and your "layer" (modification) is installing a program into it using `sudo apt install ...`. You could then use this container as the "environment" to
run particular things with the software you installed into it.
===

The Harpy team manages [a container on Dockerhub](https://hub.docker.com/repository/docker/pdimens/harpy/general) called, you guessed it, Harpy, that 
is synchronously versioned with the Harpy software. In other words, if 
you're using Harpy v1.4, it will use the container version v1.4. The 
development version of Harpy uses `latest` and the versions are automagically
managed through GitHub Actions. The Harpy container actually contains all of
the conda environments **in** it. So, when Snakemake is using the container
environment method, it will pull the versioned container from Dockerhub, and
jump in and out of container instances as required by the different jobs. 
When inside a container, Snakemake will automatically activate the correct
conda environment within the container!

## What's the Catch?
While local conda enviroments at runtime or containers might seem like  foolproof approaches, there are drawbacks.

###  Conda Caveats:
#### ⚠️ Conda Caveat 1: Inconsistent
Despite our and conda's best efforts, sometimes programs just don't install 
correctly on some systems due to unexpected system (or conda) configurations.
This results in frustrating errors where jobs fail because software that is
absolutely installed isn't being recognized (false negative), or software that wasn't 
successfully installed is being recognized (false positive). We keep trying to squash this
when it appears by using more careful dependency versioning, so improvements are being made constantly!

#### 💣 Conda Caveat 2: Troubleshooting
To manually troubleshoot many of the tasks within Harpy workflows, you
may need to jump into one of the local conda environments in `.environments`. 
That itself isn't terrible, but it's an extra step because you will
need to identify which environment is the correct one since Snakemake renames
them by their hash. There are two ways to get that information quickly:
1. check the Snakemake output log, as it lists the conda environment of the failing rule/job
2. use `harpy view environments` to list the main software in every environment available in `.environments` 

## Container Caveats
#### 🚥 Container Caveat 1: Speed
The overhead of Snakemake creating a container instance for a job, then 
cleaning it up after the job is done is not trivial and can
negatively impact runtime.

#### 💣 Container Caveat 2: Troubleshooting
The command Snakemake secretly invokes to run a job in a container is
quite lengthy. In most cases that shouldn't matter to you, but when 
something eventually goes wrong and you need to troubleshoot, it's harder
to manually rerun steps (e.g. `bwa mem genome.fa sample1.F.fq, sample1.R.fq`)
because you need a much bigger, more involved container-based command line 
call to enter a container instance and run everything with the correct
directories mounted.

