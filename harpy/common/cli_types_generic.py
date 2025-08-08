"""Module with python-click types for command-line level validations of inputs"""

import os
import yaml
import click
from pathlib import Path

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

class MultiInt(click.ParamType):
    """A class for a click type which accepts n number of integers >0 separated by commas."""
    name = "multi_int"
    def __init__(self, n: int, maximum = None, minimum = 0):
        super().__init__()
        self.n = n
        self.max = maximum
        self.min = minimum
    def convert(self, value, param, ctx):
        out = []
        parts = [i.strip() for i in value.split(',')]
        if len(parts) > self.n:
            self.fail(f"{value} is not {self.n} integers >{self.min} and separated by commas.", param, ctx)
        for i in parts:
            try:
                i = int(i)
            except ValueError:
                self.fail(f"{i} is not an integer.", param, ctx)
            if i < self.min:
                self.fail(f"All input values must be {self.min} or greater.", param, ctx)
            if self.max and i > self.max:
                self.fail(f"All input values must be {self.max} or less.", param, ctx)
            out.append(int(i))

        return out

class ReadLengths(click.ParamType):
    """A class for a click type which accepts two integers, separated by a comma."""
    name = "readlengths"
    def convert(self, value, param, ctx):
        try:
            R1,R2 = value.split(",")
        except ValueError:
            self.fail(f"{value} is not in int,int format", param, ctx)
        try:
            R1 = int(R1)
        except ValueError:
            self.fail(f"{R1} is not an integer.", param, ctx)
        try:
            R2 = int(R2)
        except ValueError:
            self.fail(f"{R2} is not an integer.", param, ctx)
        if R1 <10:
            self.fail(f"R1 reads must must be at least 10 bp. (input: {R1})", param, ctx)
        if R2 <10:
            self.fail(f"R2 reads must be at least 10 bp. (input: {R2})", param, ctx)
        return value

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

class SnakemakeParams(click.ParamType):
    """A class for a click type which accepts snakemake parameters. Does validations to make sure there isn't doubling up."""
    name = "snakemake_params"
    def convert(self, value, param, ctx):
        forbidden = "--rerun-incomplete --ri --show-failed-logs --rerun-triggers --nolock --software-deployment-method --smd --deployment --deployment-method --conda-prefix --cores -c --directory -d --snakefile -s --configfile --configfiles --conda-cleanup-pkgs --apptainer-prefix --singularity-prefix".split() 
        available = "--profile --cache --jobs -j --local-cores --resources --res --set-threads --max-threads --set-resources --set-scatter --set-resource-scopes --default-resources --default-res --preemptible-rules --preemptible-retries --envvars --touch -t --keep-going -k --force -f --executor -e --forceall -F --forcerun -R --prioritize -P --batch --until -U --omit-from -O --shadow-prefixDIR --scheduler --wms-monitor --wms-monitor-arg --scheduler-ilp-solver --conda-base-path --no-subworkflows --nosw --precommand --groups --group-components --report --report-stylesheet --reporterPLUGIN --draft-notebook --edit-notebook --notebook-listen --lint --generate-unit-tests --containerize --export-cwl --list-rules --list -l --list-target-rules --lt --dag --rulegraph --filegraph --d3dag --summary -S --detailed-summary -D --archive --cleanup-metadata --cmFILE --cleanup-shadow --skip-script-cleanup --unlock --list-changes --lc --list-input-changes --li --list-params-changes --lp --list-untracked --lu --delete-all-output --delete-temp-output --keep-incomplete --drop-metadata --version -v --printshellcmds -p --nocolor --print-compilation --force-use-threads --allow-ambiguity -a --ignore-incomplete --ii --max-inventory-time --latency-wait --output-wait -w --wait-for-files --wait-for-files-file --queue-input-wait-time --notemp --nt --all-temp --unneeded-temp-files --keep-storage-local-copies --target-files-omit-workdir-adjustment --allowed-rules --max-jobs-per-timespan --max-jobs-per-second --max-status-checks-per-second --seconds-between-status-checks --retries --restart-times -T --default-storage-provider --default-storage-prefix --local-storage-prefix --remote-job-local-storage-prefix --shared-fs-usage --scheduler-greediness --greediness --runtime-profile --local-groupid --attempt --log-handler-script --log-service --job-deploy-sources --benchmark-extended --container-image --immediate-submit --is --jobscript --js --jobname --jn --flux --container-cleanup-images --conda-not-block-search-path-envvars --conda-frontend --apptainer-args --singularity-args --use-envmodules --scheduler-solver-path --deploy-sources --target-jobs --mode --report-html-path --report-html-stylesheet-path".split()
        for i in value.split():
            if i.startswith("-"):
                if i in forbidden:
                    self.fail(f"{i} is a forbidden option because it is already used by Harpy to call Snakemake.", param, ctx)
                if i == "--dry-run" or i == "-n":
                    self.fail(f"The dry run option ({i}) is incompatible with the way Harpy launches Snakemake with a progress bar. Please use `harpy diagnose` instead.")
                if i not in available:
                    self.fail(f"{i} is not a valid Snakemake option. Run \'snakemake --help\' for a list of all Snakemake command line options.", param, ctx)
        return value

class SNPRegion(click.ParamType):
    """A class for a click type which accepts an integer, htslib-style region (chrm:start-end) or file of regions"""
    name = "snp_region"
    def convert(self, value, param, ctx):
        try:
            # is an int
            val = int(value)
            if val < 10:
                self.fail("Window size must greater than or equal to 10.", param, ctx)
            else:
                return int(value)
        except ValueError:
            pass
        if os.path.isfile(value):
            if not os.access(value, os.R_OK):
                self.fail(f"{value} is not readable. Please check file permissions and try again", param, ctx)
            with open(value, "r", encoding="utf-8") as fin:
                for idx, line in enumerate(fin, 1):
                    row = line.split()
                    if len(row) != 3:
                        self.fail(f"{value} is formatted incorrectly at line {idx}. This is the first row triggering this error, but it may not be the only one.", param, ctx)
                    else:
                        try:
                            start = int(row[1])
                            end = int(row[2])
                        except ValueError:
                            self.fail(f"{value} is formatted incorrectly at line {idx}. This is the first row triggering this error, but it may not be the only one.", param, ctx)
                    if start > end:
                        self.fail(f"The interval start position is greater than the interval end position at row {idx}. This is the first row triggering this error, but it may not be the only one.", param, ctx)
            return Path(value).resolve().as_posix()
        try:
            contig,positions = value.split(":")
        except ValueError:
            self.fail(f"{value} must be in the format contig:start-end (without spaces), where `contig` cannot contain colon (:) characters.", param, ctx)
        try:
            startpos,endpos = positions.split("-")
        except ValueError:
            self.fail(f"{value} must be in the format contig:start-end (without spaces), where `start` and `end` are integers separated by a dash (-).", param, ctx)
        try:
            startpos = int(startpos)
        except ValueError:
            self.fail(f"The region start position ({startpos}) is not a valid integer.", param, ctx)
        try:
            endpos = int(endpos)
        except ValueError:
            self.fail(f"The region end position ({endpos}) is not a valid integer.", param, ctx)
        if startpos > endpos:
            self.fail(f"The region start position ({startpos}) must be less than the end position ({endpos}).", param, ctx)
        return value

class ImputeRegion(click.ParamType):
    """A class for a click type which accepts a region in chrm:start-end-buffer format"""
    name = "snp_region"
    def convert(self, value, param, ctx):
        try:
            contig,positions = value.split(":")
        except ValueError:
            self.fail(f"{value} must be in the format contig:start-end-buffer (without spaces), where `contig` cannot contain colon (:) characters.", param, ctx)
        try:
            startpos,endpos,buffer = positions.split("-")
        except ValueError:
            self.fail(f"{value} must be in the format contig:start-end-buffer (without spaces), where `start`, `end`, and `buffer` are integers separated by a dash (-).", param, ctx)
        try:
            startpos = int(startpos)
        except ValueError:
            self.fail(f"The region start position ({startpos}) is not a valid integer.", param, ctx)
        try:
            endpos = int(endpos)
        except ValueError:
            self.fail(f"The region end position ({endpos}) is not a valid integer.", param, ctx)
        try:
            buffer = int(buffer)
        except ValueError:
            self.fail(f"The region buffer ({buffer}) is not a valid integer.", param, ctx)
        if startpos > endpos:
            self.fail(f"The region start position ({startpos}) must be less than the end position ({endpos}).", param, ctx)
        return value