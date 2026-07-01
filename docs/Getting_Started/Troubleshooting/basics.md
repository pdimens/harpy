---
label: Start Here
icon: book
order: 1
---

# :icon-book: Basic Troubleshooting
Lots of stuff can go wrong during an analysis. The intent of this page is to guide you through
navigating the inevitable errors associated with doing bioinformatics.

## Logging Streams
In shell programs, output streams are routed through standard out (`stdout`) and standard error (`stderr`)
streams. When not using `--quiet`, Harpy will output some overview text, have an updating progress bar that disappears
upon completion, and output some runtime stats at the end. Harpy will also output errors, if they occur. To make logging
Harpy consistent and pain-free, Harpy outputs the progress bar to `stdout`, while **all other text** goes to `stderr`.
What that means in practice is that if you run harpy with a `stderr` redirect to a log file, only the progress bar appears
in your console, keeping the log file clean. We do **not** recommend redirecting `stdout` (the progress bar) to a file. Here's an example:
```bash
harpy align bwa genome.fasta data/sample_*.fq.gz 2> bwa.log
```

## Troubleshooting Harpy
Harpy has two steps: first it performs checks and validations, then it runs Snakemake.

### checks and validations
First, Harpy takes your command-line inputs and checks/validates the input files and parameters.
If your parameters are not the correct type (e.g. a number where there should be a file), the
program will error immediately and inform you that something isn't right.
![Harpy option validations](/static/troubleshoot_cli.png)

If all the command-line options pass validation, then the workflow inputs
will be validated with all kinds of scripts and checks internally. For
example, an input genome FASTA will be checked to be a properly formatted
FASTA file. We do our best to include an error message that guides you
on how to fix things and try again. If there are errors you're receiving
that you think could use better wording or more information, [please let us know](https://github.com/pdimens/harpy/issues/new?assignees=&labels=enhancement&projects=&template=feature_request.yml)!
![Harpy file validations](/static/errormsg.png)

### snakemake validations
Once all the file validations pass, Harpy passes the baton over to
Snakemake. Snakemake builds a workflow graph of the rules and performs
its own checks. If you get an error before the workflow starts processing any data (there
won't yet be a Snakemake log file created), then something is wrong with
the Snakefile. Harpy may print the error to the terminal, but it's
possible you won't see any Snakemake error text (let us know!). If no
helpful error text is printed, then you should run the Snakemake command
directly and explore the output. You can copy and paste the Snakemake
command from the `config.harpy.yaml` file created by Harpy, specifically listed
under `snakemake` (either `absolute` or `relative`). If it seems like something on our end, please
[open an issue](https://github.com/pdimens/harpy/issues/new?assignees=&labels=bug&projects=&template=bug_report.yml)
and include the error text Snakemake provides.
![copy and paste this command directly to see Snakemake error text](/static/workflow_call.png)

### error during a workflow
Sometimes something goes wrong with one of the steps in a workflow. If/when
that happens, Harpy will print the offending step and all the information
Snakemake has regarding the failure. If the step had a log file, it will
print the log information too, hopefully making it easier to figure out
what's wrong. 
![error during a workflow](/static/error_text.png)
