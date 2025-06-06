#!/usr/bin/env python3

import rich_click as click
from . import align
from . import assembly, metassembly
from . import convert
from . import diagnose, resume, view
from . import deconvolve
from . import demultiplex
from . import downsample
from . import environments
from . import impute
from . import qc
from . import phase
from . import simulate
from . import snp
from . import sv
from . import template
from . import validate

click.rich_click.USE_MARKDOWN = True
click.rich_click.SHOW_ARGUMENTS = False
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = False
click.rich_click.MAX_WIDTH = 75
click.rich_click.REQUIRED_SHORT_STRING = ""
click.rich_click.ERRORS_SUGGESTION = "Try the '--help' flag for more information."
#click.rich_click.ERRORS_EPILOGUE = "Documentation: [link=https://pdimens.github.io/harpy/]https://pdimens.github.io/harpy/[/link]"

@click.group(options_metavar='', context_settings={"help_option_names" : []} )
@click.version_option("0.0.0", prog_name="harpy", hidden = True)
def cli():
    """
    An automated workflow for linked-read data
    to go from raw data to genotypes (or phased haplotypes).
    Batteries included.
    
    **demultiplex >> qc >> align >> snp >> impute >> phase >> sv**
    
    **Documentation**: [https://pdimens.github.io/harpy/](https://pdimens.github.io/harpy/)
    """

# main program
cli.add_command(align.align)
cli.add_command(assembly.assembly)
cli.add_command(convert.convert)
cli.add_command(deconvolve.deconvolve)
cli.add_command(demultiplex.demultiplex)
cli.add_command(diagnose.diagnose)
cli.add_command(downsample.downsample)
cli.add_command(environments.containerize)
cli.add_command(environments.localenv)
cli.add_command(impute.impute)
cli.add_command(metassembly.metassembly)
cli.add_command(phase.phase)
cli.add_command(qc.qc)
cli.add_command(resume.resume)
cli.add_command(simulate.simulate)
cli.add_command(snp.snp)
cli.add_command(sv.sv)
cli.add_command(validate.validate)
cli.add_command(view.view)
cli.add_command(template.template)

click.rich_click.COMMAND_GROUPS = {
    "harpy":
        [
            {
                "name": "Data Processing",
                "commands": sorted(["demultiplex","qc", "align","snp","sv","impute","phase", "simulate", "assembly", "metassembly"]),
                "panel_styles": {"border_style": "blue"}
            },
            {
                "name": "Other Commands",
                "commands": sorted(["convert","deconvolve", "downsample", "template"]),
                "panel_styles": {"border_style": "green"}
            },
            {
                "name": "Troubleshoot",
                "commands": sorted(["view", "resume", "diagnose", "validate"]),
                "panel_styles": {"border_style": "dim"}
            }
        ],
 } | align.module_docstring | convert.module_docstring | snp.module_docstring | sv.module_docstring | simulate.docstring | template.docstring

click.rich_click.OPTIONS_PANEL_TITLE = None
for i in [align, deconvolve, downsample, demultiplex, impute, phase, validate, qc, simulate, snp, sv, assembly, metassembly]:
    click.rich_click.OPTION_GROUPS |= i.docstring
