---
label: Resume
icon: arrow-right
---

When calling a workflow (e.g. [!badge corners="pill" text="qc"](qc.md)), Harpy performs various file checks
and validations, sets up the Snakemake command, output folder(s), etc. In the event you want to continue a
failed or manually terminated workflow without overwriting the workflow files (e.g. `config.harpy.yaml`),
you can use [!badge corners="pill" text="harpy resume"]. Using `resume` also skips all input/argument validations.

```bash usage
harpy resume [--conda] DIRECTORY
```

#### arguments
{.compact}
| argument    | description                                                                            |
| :---------- | :------------------------------------------------------------------------------------- |
| `DIRECTORY` | [!badge variant="info" text="required"] Output directory of an existing harpy workflow |
| `--conda`   | Generate a `/workflow/envs` folder with the necessary conda enviroments                |

The `DIRECTORY` is the output directory of a previous harpy-invoked workflow, which **must** have the `workflow/config.yaml` file.
For example, if you previously ran `harpy align bwa -o align-bwa ...`, then you would use `harpy resume align-bwa`,
which would have the necessary `workflow/config.yaml` (and other necessary things) required to successfully continue the workflow.
Using [!badge corners="pill" text="resume"] does **not** overwrite any preprocessing files in the target directory (whereas rerunning the workflow would),
which means you can also manually modify the `config.yaml` file (advanced, not recommended unless you are confident with what you're doing).

[!badge corners="pill" text="resume"] also requires an existing and populated `workdir/envs/` directory in the target directory, like the kind all
main `harpy` workflows would create. If one is not present, you can use `--conda` to create one.
