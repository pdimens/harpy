#!/usr/bin/env python3

import rich_click as click
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
from . import resume
from . import deconvolve
from .popgroups import popgroup
from .stitchparams import stitchparams

click.rich_click.USE_MARKDOWN = True
click.rich_click.SHOW_ARGUMENTS = False
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = False
click.rich_click.MAX_WIDTH = 75
click.rich_click.REQUIRED_SHORT_STRING = ""
click.rich_click.ERRORS_SUGGESTION = "Try the '--help' flag for more information."
click.rich_click.ERRORS_EPILOGUE = """
See the documentation: [link=https://pdimens.github.io/harpy/]https://pdimens.github.io/harpy/[/link]
"""

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
@click.version_option("0.0.0", prog_name="harpy")
def cli():
    """
    An automated workflow for haplotagging linked-read data
    to go from raw data to genotypes (or phased haplotypes).
    Batteries included.
    
    **demultiplex >> qc >> align >> snp >> impute >> phase >> sv**
    
    **Documentation**: [https://pdimens.github.io/harpy/](https://pdimens.github.io/harpy/)
    """

# main program
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
cli.add_command(resume.resume)
cli.add_command(deconvolve.deconvolve)

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
                "commands": ["resume", "hpc", "preflight", "deconvolve", "popgroup", "stitchparams"]
            }
        ],
 } | simulate.commandstring | hpc.docstring

for i in [align, deconvolve, demultiplex, impute, phase, preflight, qc, simulate, snp, sv]:
    click.rich_click.OPTION_GROUPS |= i.docstring