---
label: View
icon: eye
---

This convenience command lets you view relevant details within a Harpy workflow directory
without having to fish around for the right files.

```bash usage
harpy view MODE DIRECTORY
```

```bash example | view the workflow configuration of an existing harpy-generated folder
harpy view config Align/bwa
```

### arguments
{.compact}
| argument      |  description                                                                           |
| :------------ | :------------------------------------------------------------------------------------- |
| `DIRECTORY`   | [!badge variant="info" text="required"] Output directory of an existing harpy workflow |

### modes
{.compact}
| MODE            | description                                                       |
|:----------------|:------------------------------------------------------------------|
| `config`        | View a workflow's config file                                     |
| `environments`* | View the conda environments and their software in `.environments` |
| `log`           | View a workflow's last log file                                   |
| `snakefile`     | View a workflow's snakefile                                       |
| `snakeparams`   | View a workflow's snakemake parameter file                        |

#### environments
`view environments` is an exception in that it does not require a `DIRECTORY` argument.
To use `environments`, you can run it without arguments or give it a `SOFTWARE` argument
to only print environments where `SOFTWARE` is found (it also works with partial matches).

``` viewing all environments
harpy view environments

.environments/f5cb053d77e72fca1c7b6463448fd855_
  - falco=1.2.4
  - fastp
  - multiqc=1.28
  - pysam=0.22
  - quickdeconvolution

.environments/9f9995de2cc77c81654693b8fb002922_
  - quarto
  - r-dt
  - r-dplyr
  - r-highcharter
  - r-magrittr
  - r-plotly
  - r-scales
  - r-tidyr
  - r-viridislite
  - r-xml2
  - r-biocircos
```

``` searching for a program
harpy view environments fast

.environments/f5cb053d77e72fca1c7b6463448fd855_
  - falco=1.2.4
  → fastp
  - multiqc=1.28
  - pysam=0.22
  - quickdeconvolution
```