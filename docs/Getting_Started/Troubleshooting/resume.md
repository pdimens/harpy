---
label: Resume
icon: arrow-right
---

# :icon-arrow-right: Resume

When calling a workflow (e.g. [!badge corners="pill" text="qc"](/Commands/qc.md)), Harpy performs various file checks
and validations, sets up the Snakemake command, output folder(s), etc. In the event you want to continue a
failed or manually terminated workflow without overwriting the workflow files (e.g. `workflow.yaml`),
you can use [!badge corners="pill" text="harpy resume"]. Using `resume` also bypasses all input/argument validations.

```bash usage
harpy resume [--threads] [--direct] DIRECTORY
```

#### arguments
{.compact .clean}
| argument {.whitespace-nowrap} | description                                                                                                               |
| :---------------------------- | :------------------------------------------------------------------------------------------------------------------------ |
| `DIRECTORY`                   | [!badge variant="info" text="required"] Output directory of an existing harpy workflow                                    |
| `--conda`                     | [!badge variant="danger" text="being deprecated"] Generate a `/workflow/envs` folder with the necessary conda enviroments |
| `--direct`                    | [!badge text="unreleased"] Launch Snakemake without any Harpy intervention                                                |
| `--threads`                   | Change the number of threads the workflow will be run with                                                                |

The `DIRECTORY` is the output directory of a previous harpy-invoked workflow, which **must** have the `workflow/config.yaml`
and `workflow/workflow.yaml` files.
For example, if you previously ran `harpy align bwa -o align-bwa ...`, then you would use `harpy resume align-bwa`,
which would have the necessary `workflow/config.yaml` (and other necessary things) required to successfully continue the workflow.
Using [!badge corners="pill" text="resume"] does **not** overwrite any preprocessing files in the target directory (whereas rerunning the workflow would),
which means you can also manually modify the `config.yaml` file (advanced, not recommended unless you are confident with what you're doing).

[!badge corners="pill" text="resume"] also requires an existing and populated `workdir/envs/` directory in the target directory, like the kind all
main `harpy` workflows would create. If one is not present, you can use `--conda` to create one (being deprecated).

## Considerations
The snakefiles in harpy workflows are, by design, not strict for the presence/absence of `Parameter` keys in a workflow's
corresponding `workflow.yaml` file. If a parameter key is absent, the workflow will default to using that parameter's
CLI default value. This silent behavior can be considered both a bug and a feature. For example, here is the expected
parameter section of `workflow.yaml` from `harpy preprocess meier2021`:
```yaml
Workflow:
    ...
Parameters:
  qx-rx: true
  unknown-barcodes: false
  unknown-samples: true
  stitch:
    base: false
    complementary: false
Inputs:
    ...
```
The command line defaults for each of these parameters is `false`, meaning a parameter section like this would result
in the same workflow:
```yaml
Workflow:
    ...
Parameters:
  qx-rx: true
  unknown-samples: true
Inputs:
    ...
```

In most cases, workflows are started using standard Harpy commands like `harpy align bwa`, which guarantees correct
`workflow.yaml` files ingested by Snakemake. The two typical use-cases of `resume` are to restart a workflow that cutoff
midway or to initiate a modified workflow without validations, and we can't guard against one behavior without
spamming notices or errors for the other.  We decided that hand-editing `workflow.yaml` files (or snakefiles) will be
considered an advanced **at your own risk** use-case, and Harpy will not inform you that parameters keys are missing so as to not
be an obstacle to customizing workflows.