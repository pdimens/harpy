---
label: Development
icon: tools
order: 1
---

# :icon-tools: Developing Harpy
Harpy is an open source program written using a combination of BASH, R, 
RMarkdown, Python, and Snakemake. This page provides information on Harpy's
development and how to contribute to it, if you were inclined to do so.

=== Design Philosophy

Before we get into the technical details, you, dear reader, need to understand
why Harpy is the way it is. Harpy may be a pipeline for other software, but 
there is a lot of extra stuff built in to make it user 
friendly. Not just friendly, but _compassionate_. The guiding ethos for Harpy is
**"We don't hate the user"**. That means there is a lot
of code that checks input files, runtime details, etc. to exit before 
Snakemake takes over. This is done to minimize time wasted on minor 
errors that only show their ugly heads 18 hours into a 96 hour process. With that in mind:
1. **Code should be written clearly** because someone else will need to read it at 
some point, and that person could be future-you who hasn't seen or thought 
about that code for a while. Write nicely. Annotate.
2. **Error messages should provide all the information a user needs to fix the problem and retry**. It's not enough to exit when an error is identified. Collate
the things causing the error, explain to the user what and why. Harpy follows a
style of presenting and explaining the error, then providing a solution and showing exactly what files/rows/columns/etc. caused the error. Be kind to users.

![These are Harpy error messages](/static/errormsg.png)

3. **Documentation is just as important as the code**. No features are undocumented,
and the documentation should read like something that a new student can
pick up and understand. Good documentation, written compassionately, will lower
the barrier of entry to people who just want to process their haplotag data. Harpy
isn't about ego, it's about accessibility. We invest in this documentation because
we _want_ to, not because we need to.
===

## Installing Harpy for development
The process follows cloning the harpy repository, installing the preconfigured conda environment, and running the `resources/buildlocal.sh`
script to move all the necessary files to the `/bin/` path within your active conda environment.

==- Step 1: clone the repository

```bash clone the repository
git clone https://github.com/pdimens/harpy.git
```

==- Step 2: install the conda environment dependencies
```bash install the dependencies with conda/mamba
mamba env create --name harpy --file resources/harpy.yaml
```
This will create a conda environment named `harpy` with all the bits necessary to successfully run Harpy. You can change the name of this environment by specifying
`--name something`. 

==- Step 3: activate the environment
The environment with all the preinstalled dependencies can be activated with:
```bash activate the conda environment
# assuming the environment name is harpy from the step above
mamba activate harpy
```
==- Step 4: install the Harpy files

Call the `resources/buildlocal.sh` bash script to finish the installation.
This will build the `harpy` python program, and copy all the additional files Harpy needs to run
to the `bin/` directory of your active conda environment.
```bash install harpy and the necessary files
bash resources/buildlocal.sh
```
===

## Harpy's components
### source code
Harpy runs in two stages:
1. it recieves command line inputs and parses them
2. uses the parsed command line inputs to run a specific Snakemake workflow

To accomplish this, Harpy is written as a Python program, using [rich-click](https://github.com/ewels/rich-click)
for the aesthetically pleasing interface. Since Harpy also relies on Snakemake
snakefiles for each module, bash scripts, python scripts, and rmarkdown files,
not all of it can be installed as a pure python program using `setuptools`.
The build process installs part of Harpy as a pure-python command line program, but
all the extra files Harpy needs to run need to be installed separately. All of 
this is handled by `resources/buildlocal.sh`. It's a little circuitous, but it's how
we can keep the source code modular, installable, and have the flexibility of 
using non-python code.

### Bioconda recipe
For the ease of installation for end-users, Harpy has a [recipe](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/harpy/meta.yaml) 
and [build script](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/harpy/build.sh) in Bioconda,
which makes it available for download and installation. A copy of the recipe and
build script is also stored in `resources/meta.yml` and `resources/build.sh`. The yaml file
is the metadata of the package, including software deps and their versions. The
build script is how conda will install all of Harpy's parts. Now that Harpy is hosted on
bioconda, when a new version is tagged with a release, Bioconda will automatically create
a pull request (after a delay), typically not requiring any intervention on the development side
for the newest Harpy version to be released for conda installation.


## The Harpy repository
### structure
Harpy exists as a Git repository and has 5 standard branches that are used in 
specific ways during development. Git is a popular version control system and 
discussing its use is out of the scope of this documentation, however there is no 
shortage of [great resources](https://www.youtube.com/watch?v=8Dd7KRpKeaE) to get 
you started. The 5 standard branches in the Harpy repository are outlined in the 
table below:

| branch | purpose |
|:---| :---|
| `main` | staging and testing area for new code prior to creating the next release  |
| `docs` | the source documentation files (markdown and configs) that are deployed for the current documentation |
| `docs_dev` | the updated documentation intended for the next release and is not deployed |
| `retype` | the branch that `docs` deploys to and contains the current rendered documentation and is not to be touched |

### development workflow
The dev workflow is reasonably standard:
1. [create a fork](https://github.com/pdimens/harpy/fork) of Harpy
2. within your fork, create a new branch, name it something relevant to what you intend to do (_e.g._, `naibr_bugfix`, `add_deepvariant`)
3. add and modify code with your typical coding workflow, pushing your changes to your Harpy fork
4. when it's ready for inclusion into Harpy (and testing), create a Pull Request to merge your changes into the Harpy `main` branch

### containerization
As of Harpy v1.0, the software dependencies that the Snakemake workflows use are pre-configured as a Docker image
that is uploaded to Dockerhub. Updating or editing this container can be done automatically or manually.

#### automatically
The [rebuildcontainer.yml](https://github.com/pdimens/harpy/blob/main/.github/workflows/rebuildcontainer.yml)
GitHub Action will automatically create a Dockerfile with [!badge corners="pill" text="harpy containerize"] (a hidden harpy command)
 and build a new Docker container, then upload it to [dockerhub](https://hub.docker.com/repository/docker/pdimens/harpy/general)
with the `latest` tag. This action is triggered on `push` or `pull request` with changes to either
`src/harpy/conda_deps` or `src/harpy/snakefiles/containerize.smk`.

#### manually
The dockerfile for that container is created by using a hidden harpy command [!badge corners="pill" text="harpy containerize"]
```bash auto-generate Dockerfile
harpy containerize
```
which does all of the work for us. The result is a `Dockerfile` that has all of the conda environments
written into it. After creating the `Dockerfile`, the image must then be built.

```bash build the Docker image
cd resources
docker build -t pdimens/harpy .
```
This will take a bit because the R dependencies are hefty. Once that's done, the image can be pushed to Dockerhub:
```bash push image to Dockerhub
docker push pdimens/harpy
```
This containerize -> dockerfile -> build -> process will push the changes to Dockerhub with the `latest` tag, which is suitable for
the development cycle. When the container needs to be tagged to be associated with the release of a new Harpy version, you will need to
add a tag to the `docker build` step:
```bash build tagged Docker image
cd resources
docker build -t pdimens/harpy:TAG
```
where `TAG` is the Harpy version, such as `1.0`, `1.4.1`, `2.1`, etc. As such, during development, the `containerized: docker://pdimens/harpy:TAG` declaration at the top of the snakefiles should use the `latest` tag, and when ready for release, changed to match the Harpy
version. So, if the Harpy version is `1.4.12`, then the associated docker image should also be tagged with `1.4.12`. The tag should remain `latest`
(unless there is a very good reason otherwise) since automatic Docker tagging happens upon releases of new Harpy versions.

## Automations
### Testing
CI (**C**ontinuous **I**ntegration) is a term describing automated actions that do
things to/with your code and are triggered by how you interact with a repository.
Harpy has a series of GitHub Actions triggered by interactions with the `main` branch (in `.github/workflows`) 
to test the Harpy modules depending on which files are being changed by the push or
pull request. It's setup such that, for example, when files associated with 
demultiplexing are altered, it will run `harpy demultiplex` on the test data 
in the cloud and notify the Harpy devs if for some reason `harpy demultiplex`
could not run successfully to completion. These tests do not test for accuracy,
but test for breaking behavior. You'd be shocked to find out how many errors
crop up this way and require more development so Harpy can be resilient to more use cases.
### Releases
There is [an automation](https://github.com/pdimens/harpy/blob/dev/.github/workflows/createrelease.yml)
that gets triggered every time Harpy is tagged with the new version. It strips out the unnecessary files and will
upload a cleaned tarball to the new release (reducing filesize by orders of magnitude). The automation will also
build a new Dockerfile and tag it with the same git tag for Harpy's next release and push it to Dockerhub.
In doing so, it will also replace the tag of the container in all of Harpy's snakefiles from `latest` to the
current Harpy version. In other words, during development the top of every snakefile reads
`containerized: docker://pdimens/harpy:latest` and the automation replaces it with (e.g.) `containerized: docker://pdimens/harpy:1.17`.
Tagging is easily accomplished with Git commands in the command line:
```bash
# make sure you're on the main branch
$ git checkout main

# create the tag locally, where X.X.X is something like 1.7.1
$ git tag X.X.X

# push the new tag to the repository
$ git push origin X.X.X
```
