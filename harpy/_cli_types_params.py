"""Module with python-click types for command-line level validations of program-specific extra-params inputs"""

import os
import click

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

