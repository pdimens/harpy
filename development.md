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
The process follows cloning the harpy repository, installing the preconfigured conda environment, and running the `misc/buildlocal.sh`
script to move all the necessary files to the `/bin/` path within your active conda environment.

==- Step 1: clone the repository

```bash clone the repository
git clone https://github.com/pdimens/harpy.git
```

==- Step 2: install the conda environment dependencies
```bash install the dependencies with conda/mamba
mamba env create --name harpy --file misc/harpyenv.yaml
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

Call the `misc/buildlocal.sh` bash script to finish the installation.
This will build the `harpy` python program, and copy all the additional files Harpy needs to run
to the `bin/` directory of your active conda environment.
```bash install harpy and the necessary files
bash misc/buildlocal.sh
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
this is handled by `misc/buildlocal.sh`. It's a little circuitous, but it's how
we can keep the source code modular, installable, and have the flexibility of 
using non-python code.

### Bioconda recipe
For the ease of installation for end-users, Harpy has a [recipe](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/harpy/meta.yaml) 
and [build script](https://github.com/bioconda/bioconda-recipes/blob/master/recipes/harpy/build.sh) in Bioconda,
which makes it available for download and installation. A copy of the recipe and
build script is also stored in `misc/meta.yml` and `misc/build.sh`. The yaml file
is the metadata of the package, including software deps and their versions. The
build script is how conda will install all of Harpy's parts. In order to modify 
these files for a new release, you need to fork `bioconda/bioconda-recipes`, 
create a new branch, modify the Harpy `meta.yml` (and possibly `build.sh`) files. Bioconda
has an bot that looks for changes to the version number in the `meta.yml` file
and will automatically submit a Pull Request when it notices that's been changed.

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
| `main` | the source code of the current release and used for new bioconda releases |
| `dev`  | staging and testing area for new code prior to merging with `main` for release |
| `docs` | the source documentation files (markdown and configs) that are deployed for the current documentation |
| `docs_dev` | the updated documentation intended for the next release and is not deployed |
| `retype` | the branch that `docs` deploys to and contains the current rendered documentation and is not to be touched |

### development workflow
The dev workflow is reasonably standard:
1. [create a fork](https://github.com/pdimens/harpy/fork) of Harpy
2. within your fork, create a new branch, name it something relevant to what you intend to do (_e.g._, `naibr_bugfix`, `add_deepvariant`)
    - create a branch off of `main` if you are trying to fix a bug in the release version
    - create a branch off of `dev` if you are adding a new feature or breaking change
3. add and modify code with your typical coding workflow, pushing your changes to your Harpy fork
4. when it's ready for inclusion into Harpy (and testing), create a Pull Request to merge your changes into the Harpy `dev` branch

## Testing and CI
CI (**C**ontinuous **I**ntegration) is a term describing automated actions that do
things to/with your code and are triggered by how you interact with a repository.
Harpy has a series of GitHub Actions triggered by interactions with the `dev` branch (in `.github/workflows`) 
to test the Harpy modules depending on which files are being changed by the push or
pull request. It's setup such that, for example, when files associated with 
demultiplexing are altered, it will run `harpy demultiplex` on the test data 
in the cloud and notify the Harpy devs if for some reason `harpy demultiplex`
could not run successfully to completion. These tests do not test for accuracy,
but test for breaking behavior. You'd be shocked to find out how many errors
crop up this way and require more development so Harpy can be resilient to more use cases.