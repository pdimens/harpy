"""Module with python-click types for command-line level validations of inputs"""

import os
import click

class IntList(click.ParamType):
    """A class for a click type which accepts an arbitrary number of integers, separated by a comma."""
    name = "int_list"
    def __init__(self, entries):
        super().__init__()
        self.entries = entries

    def convert(self, value, param, ctx):
        try:
            parts = [i.strip() for i in value.split(',')]
            if len(parts) != self.entries:
                raise ValueError
            for i in parts:
                try:
                    int(i)
                except:
                    raise ValueError
            return [int(i) for i in parts]
        except ValueError:
            self.fail(f"{value} is not a valid list of integers. The value should be {self.entries} integers separated by a comma.", param, ctx)

class KParam(click.ParamType):
    """A class for a click type which accepts any number of odd integers separated by a comma, or the word auto."""
    name = "k_param"
    def convert(self, value, param, ctx):
        try:
            if value == "auto":
                return value
            parts = [i.strip() for i in value.split(',')]
            for i in parts:
                if int(i) % 2 == 0 or int(i) > 128:
                    raise ValueError
            # make it a set as a failsafe against duplicates
            return list(set(int(i) for i in parts))
        except ValueError:
            self.fail(f"{value} is not 'auto' or odd integers <128 separated by a comma.", param, ctx)

class ContigList(click.ParamType):
    """A class for a click type which accepts a file of contigs or a list of contigs separated by a comma."""
    name = "contig_list"
    def convert(self, value, param, ctx):
        # check if it's a file
        if os.path.exists(value):
            if not os.path.isfile(value):
                self.fail(f"{value} is not a file.", param, ctx)
            with open(value, "r") as cont_in:
                return [i.strip() for i in cont_in.readlines()]
        else:
            # make it a set as a failsafe against duplicates
            return list(set(i.strip() for i in value.split(',')))

class InputFile(click.ParamType):
    """A class for a click type that verifies that a file exists and that it has an expected extension"""
    name = "input_file"
    def __init__(self, filetype, gzip_ok):
        super().__init__()
        self.filetype = filetype
        self.gzip_ok = gzip_ok
    def convert(self, value, param, ctx):
        filedict = {
            "fasta": [".fasta", ".fas", ".fa", ".fna", ".ffn", ".faa", ".frn"],
            "vcf": ["vcf", "bcf", "vcf.gz"],
            "gff": [".gff",".gff3"]
        }
        if self.filetype not in filedict:
            self.fail(f"Extension validation for {self.filetype} is not yet implemented. This error should only appear during development; if you are a user and seeing this, please post an issue on GitHub: https://github.com/pdimens/harpy/issues/new?assignees=&labels=bug&projects=&template=bug_report.yml")
        if not os.path.exists(value):
            self.fail(f"{value} does not exist. Please check the spelling and try again.", param, ctx)
        elif not os.access(value, os.R_OK):
            self.fail(f"{value} is not readable. Please check file/directory permissions and try again", param, ctx)
        if os.path.isdir(value):
            self.fail(f"{value} is a directory, but input should be a file.", param, ctx)
        valid = False
        lowercase = value.lower()
        for ext in filedict[self.filetype]:
            valid = True if lowercase.endswith(ext) else valid
            if self.gzip_ok:
                valid = True if lowercase.endswith(ext + ".gz") else valid
        if not valid and not self.gzip_ok:
                self.fail(f"{value} does not end with one of the expected extensions [" + ", ".join(filedict[self.filetype]) + "]. Please verify this is the correct file type and rename the extension for compatibility.", param, ctx)
        if not valid and self.gzip_ok:
            self.fail(f"{value} does not end with one of the expected extensions [" + ", ".join(filedict[self.filetype]) + "]. Please verify this is the correct file type and rename the extension for compatibility. Gzip compression (ending in .gz) is allowed.", param, ctx)
        return value

class SnakemakeParams(click.ParamType):
    """A class for a click type which accepts snakemake parameters. Does validations to make sure there isn't doubling up."""
    name = "snakemake_params"
    def convert(self, value, param, ctx):
        forbidden = "--rerun-incomplete --ri --show-failed-logs --rerun-triggers --nolock --software-deployment-method --smd --deployment --deployment-method --conda-prefix --cores -c --directory -d --snakefile -s --configfile --configfiles".split() 
        available = "--dry-run --dryrun -n --profile --cache --jobs -j --local-cores --resources --res --set-threads --max-threads --set-resources --set-scatter --set-resource-scopes --default-resources --default-res --preemptible-rules --preemptible-retries --envvars --touch -t --keep-going -k --force -f --executor -e --forceall -F --forcerun -R --prioritize -P --batch --until -U --omit-from -O --shadow-prefixDIR --scheduler --wms-monitor --wms-monitor-arg --scheduler-ilp-solver --conda-base-path --no-subworkflows --nosw --precommand --groups --group-components --report --report-stylesheet --reporterPLUGIN --draft-notebook --edit-notebook --notebook-listen --lint --generate-unit-tests --containerize --export-cwl --list-rules --list -l --list-target-rules --lt --dag --rulegraph --filegraph --d3dag --summary -S --detailed-summary -D --archive --cleanup-metadata --cmFILE --cleanup-shadow --skip-script-cleanup --unlock --list-changes --lc --list-input-changes --li --list-params-changes --lp --list-untracked --lu --delete-all-output --delete-temp-output --keep-incomplete --drop-metadata --version -v --printshellcmds -p --debug-dag --nocolor --quiet -q --print-compilation --verbose --force-use-threads --allow-ambiguity -a --ignore-incomplete --ii --max-inventory-time --latency-wait --output-wait -w --wait-for-files --wait-for-files-file --queue-input-wait-time --notemp --nt --all-temp --unneeded-temp-files --keep-storage-local-copies --target-files-omit-workdir-adjustment --allowed-rules --max-jobs-per-timespan --max-jobs-per-second --max-status-checks-per-second --seconds-between-status-checks --retries --restart-times -T --wrapper-prefix --default-storage-provider --default-storage-prefix --local-storage-prefix --remote-job-local-storage-prefix --shared-fs-usage --scheduler-greediness --greediness --no-hooks --debug --runtime-profile --local-groupid --attempt --log-handler-script --log-service --job-deploy-sources --benchmark-extended --container-image --immediate-submit --is --jobscript --js --jobname --jn --flux --container-cleanup-images --use-conda --conda-not-block-search-path-envvars --list-conda-envs --conda-cleanup-envs --conda-cleanup-pkgs --conda-create-envs-only --conda-frontend --use-apptainer --use-singularity --apptainer-prefix --singularity-prefix --apptainer-args --singularity-args --use-envmodules --scheduler-solver-path --deploy-sources --target-jobs --mode --report-html-path --report-html-stylesheet-path".split()
        for i in value.split():
            if i.startswith("-"):
                if i in forbidden:
                    self.fail(f"{i} is a forbidden option because it is already used by Harpy to call Snakemake.", param, ctx)
                if i not in available:
                    self.fail(f"{i} is not a valid Snakemake option. Run \'snakemake --help\' for a list of all Snakemake command line options.", param, ctx)
        return value

class HPCProfile(click.ParamType):
    """A class for a click type which accepts a directory with a snakemake HPC profile. Does validations to make sure the config file is there."""
    name = "hpc_profile"
    def convert(self, value, param, ctx):
        if not os.path.exists(value):
            self.fail(f"{value} does not exist. Please check the spelling and try again.", param, ctx)
        elif not os.access(value, os.R_OK):
            self.fail(f"{value} is not readable. Please check file/directory permissions and try again", param, ctx)
        if os.path.isfile(value):
            self.fail(f"{value} is a file, but input should be a directory.", param, ctx)
        if not os.path.exists(f"{value}/config.yaml"):
            self.fail(f"{value} does not contain the necessary config.yaml file.", param, ctx)
        elif not os.access(f"{value}/config.yaml", os.R_OK):
            self.fail(f"{value}/config.yaml does not have read access. Please check the file permissions and try again.", param, ctx)
        return value
    