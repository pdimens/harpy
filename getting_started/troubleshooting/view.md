# :icon-eye: View

This convenience command lets you view relevant details within a Harpy workflow directory
without having to fish around for the right files. Since version 4.0, a harpy installation
also aliases `harpy view` as `hv` for convenience.

```bash usage
harpy view MODE DIRECTORY
# or #
hv MODE DIRECTORY
```

```bash example | view the workflow configuration of an existing harpy-generated folder
harpy view config Align/bwa
```

#### The `.harpyerror` file
A feature introduced in v4 is the automatic creation of a `.harpyerror` file in the work
directory. This file contains a string path to the last output directory that caused an error.
When any of the `view` commands are run without the `DIRECTORY` argument (except `envs`),
harpy will look for the `.harpyerror` 

### arguments
{.compact .clean}
| argument    | description                                                                                              |
| :---------- | :------------------------------------------------------------------------------------------------------- |
| `DIRECTORY` | Output directory of an existing harpy workflow. Defaults to path in `.harpyerror` if present and omitted |

### modes
{.compact .clean}
| Mode   {.whitespace-nowrap} | description                                                       |
| :-------------------------- | :---------------------------------------------------------------- |
| `config`                    | View a workflow's config file                                     |
| `error`                     | View the exact failture of a snakemake run                        |
| `envs`*                     | View the conda environments and their software in `.environments` |
| `log`                       | View a workflow's last log file                                   |
| `snakefile`                 | View a workflow's snakefile                                       |
| `profile`                   | View a workflow's snakemake parameter file                        |

#### environments
`view envs` is an exception in that it does not require a `DIRECTORY` argument.
To use `environments`, you can run it without arguments or give it a `SOFTWARE` argument
to only print environments where `SOFTWARE` is found (it also works with partial matches).

``` viewing all environments
harpy view environments

Conda Environments
├── .environments/e38c41849c266f1765a6dad2ea0db37a_
│   ├── bwa
│   ├── minibwa
│   ├── minimap2
│   ├── samtools=1.23
│   ├── seqtk
│   ├── strobealign
│   └── tabix
├── .environments/4dc2b26a762f0f46c3c195218651080a_
│   ├── dmox>=0.2
│   └── pheniqs=2.1
└── .environments/0ffc3ee8f456a1920afe71babd002bf7_
    ├── bcftools=1.23
    ├── freebayes=1.3.9
    ├── leviathan
    └── naibr-plus=0.5.4
```

``` searching for a program
harpy view environments fast

Conda Environments
└── .environments/3090e313a5f362f817f7572c06ff912a_
    ├── click=8.2.1
    ├── falco=1.2.5
    ├── fastp <
    ├── mosdepth
    ├── multiqc=1.30
    ├── pysam=0.23
    ├── quickdeconvolution
    └── samtools
```
