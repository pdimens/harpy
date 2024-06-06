#!/usr/bin/env python3

import os
import sys

try:
    harpypath = '{CONDA_PREFIX}'.format(**os.environ) + "/bin"
except:
    enverror = "\033[1;33mERROR:\033[00m Harpy expects to run from within an active conda environment, but one was not detected."
    print(enverror, file = sys.stderr)
    fix = "\033[1;34mSOLUTION:\033[00m Activate the conda environment Harpy was installed into and run Harpy again."
    print("\n" + fix, file = sys.stderr)
    print(f"\n\033[1mDetails:\033[00m")
    details = "In order to work correctly, Harpy expects several software packages to be available in the PATH, which are provided automatically with Harpy's conda-based installation. If you're seeing this message, no active conda environment was detected upon execution, and Harpy exited as an early failsafe against unexpected runtime errors associated with \"missing\" files and packages."
    print(details, file = sys.stderr)
    sys.exit(1)

from . import align
from . import demultiplex
from . import impute
from . import phase
from . import preflight
from . import qc
from . import simulate
from . import snp
from . import sv
from . import container
from . import hpc
from .popgroups import popgroup
from .stitchparams import stitchparams
import rich_click as click

click.rich_click.USE_MARKDOWN = True
click.rich_click.SHOW_ARGUMENTS = False
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = False
click.rich_click.MAX_WIDTH = 75
click.rich_click.REQUIRED_SHORT_STRING = ""
click.rich_click.ERRORS_SUGGESTION = "Try the '--help' flag for more information."
click.rich_click.ERRORS_EPILOGUE = "See the documentation: [link=https://pdimens.github.io/harpy/]https://pdimens.github.io/harpy/[/link]"

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
@click.version_option("1.1.0", prog_name="Harpy")
def cli():
    """
    ## Harpy haplotagging pipeline
    
    An automated workflow to demultiplex sequences, trim and qc reads, 
    map sequences, call variants, impute genotypes, and phase 
    haplotypes of Haplotagging data. Batteries included.
    
    **demultiplex >> qc >> align >> snp >> impute >> phase >> sv**
    
    **Documentation**: [https://pdimens.github.io/harpy/](https://pdimens.github.io/harpy/)
    """

# main program
#cli.add_command(hpc)
cli.add_command(popgroup)
cli.add_command(stitchparams)
cli.add_command(preflight.preflight)
cli.add_command(demultiplex.demultiplex)
cli.add_command(qc.qc)
cli.add_command(align.align)
cli.add_command(snp.snp)
cli.add_command(sv.sv)
cli.add_command(impute.impute)
cli.add_command(phase.phase)
cli.add_command(simulate.simulate)
cli.add_command(container.containerize)
cli.add_command(hpc.hpc)

## the modules ##
click.rich_click.COMMAND_GROUPS = {
    "harpy":
        [
            {
                "name": "Modules",
                "commands": ["demultiplex","qc", "align","snp","sv","impute","phase", "simulate"],
            },
            {
                "name": "Other Commands",
                "commands": ["hpc", "preflight", "popgroup", "stitchparams"]
            }
        ],
 } | simulate.commandstring | hpc.docstring

click.rich_click.OPTION_GROUPS = demultiplex.docstring | preflight.docstring | qc.docstring | align.docstring | snp.docstring | sv.docstring | impute.docstring | phase.docstring | simulate.docstring
