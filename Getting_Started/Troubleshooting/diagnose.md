---
label: Diagnose
icon: unverified
---

# :icon-unverified: harpy diagnose

Harpy includes a series of [!badge corners="pill" text="diagnose"] commands to faciliate troubleshooting


## diagnose stall
The [!badge corners="pill" text="diagnose stall"] command will run dry-run Snakemake for a workflow with the `--dry-run --debug-dag` options
and print output to try to identify why Snakemake might be stalling for a workflow.
This is typically used during development or for troubleshooting edge-cases and you're likely/hopefully
not going to need to use this feature. Similar to [!badge corners="pill" text="resume"](/Getting_Started/Troubleshooting/resume.md),
this requires a Harpy-generated output folder as the sole argument.

```bash usage
harpy diagnose stall path/to/dir
```

```bash example
harpy diagnose stall Align/bwa
```

## diagnose rule [!badge varaint="warning" text="unreleased"]
The [!badge corners="pill" text="diagnose rule"] command is a shortcut for trying to manually rerun a failing rule.The
command will scan the input directory for the most recent Snakemake log file and parse it to identify the first rule
triggering the workflow to fail. The input files of the failing rule will be checked for their existence, and if not,
Harpy will call Snakemake directly with `--no-temp` and the missing files as targets to ensure the creation of the missing
input files (this is skipped if all input files are present). After, the failing command will be run directly, without invoking Snakemake, in any appropriate container or conda environment.

```bash usage
harpy diagnose rule path/to/dir
```

```bash example
harpy diagnose rule Deconvolve/
```

The output walks you through what Harpy is trying to do:
![terminal output of diagnose rule](/static/diagnose_rule.png)