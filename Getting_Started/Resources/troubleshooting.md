---
label: Troubleshooting
icon: alert
order: 5
---

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

---
## Common Issues
### installation issue
Conda is an awesome package manager, but _was_ slow and used a ton of memory
as dependencies increased. Recent (`23.10+`) versions of Conda [now use the `libmamba` solver](https://www.anaconda.com/blog/a-faster-conda-for-a-growing-community),
the super-fast and super-lightweight solver from Mamba. If you're experiencing
a suspiciously slow Harpy installation, either update Conda to at least version `23.10` or use Mamba.

### imputation or phasing failure
If you use `bamutils clipOverlap` on alignments that are used for the [!badge corners="pill" text="impute"](/Workflows/impute.md) or
[!badge corners="pill" text="phase"](/Workflows/phase.md) modules, they will cause both programs to error. We don't know why, but they do.

**Solution**: Do not clip overlapping alignments for bam files you intend to use for
the [!badge corners="pill" text="impute"](/Workflows/impute.md) or
[!badge corners="pill" text="phase"](/Workflows/phase.md) modules. Harpy does not clip overlapping alignments, so
alignments produced by Harpy should work just fine.

### SAM name and ID mismatch
Aligning a sample to a genome via Harpy will insert the sample name (based on the file name)
into the alignment header (the `@RG ID:name SM:name` tag). It likewise expects, through various steps,
that the sample names in resulting vcf files match the filenames of associated bam files. This creates 
problems when manually renaming alignment files after the creation of any vcf files. If you rename the 
bam file, the alignments will still have the original sample name hardcoded into the file header. 
Harpy will check for this and will preemtively warn you of a mismatch between file name and encoded
sample name. Due to certain expectations of the workflow, this mismatch will absolutely cause things
to fail, hence the pre-flight check.

**Solution**: If you need to rename a bam file, do so using the [rename_bam](utilities.md#rename_bam) script bundled with Harpy, which is a just a thin veneer over `samtools addreplacerg` with some extra validations.
```bash
rename_bam newname input.bam 
```
Call the script with no arguments to see the full usage instructions.


**More cases will be added here as they become apparent to us**
