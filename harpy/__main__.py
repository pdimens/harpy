#!/usr/bin/env python3

import rich_click as click
from harpy.commands import align
from harpy.commands import assembly, metassembly
from harpy.commands import convert
from harpy.commands import diagnose, resume, view
from harpy.commands import deconvolve
from harpy.commands import demultiplex
from harpy.commands import downsample
from harpy.commands import environments
from harpy.commands import impute
from harpy.commands import qc
from harpy.commands import phase
from harpy.commands import simulate
from harpy.commands import snp, sv
from harpy.commands import template
from harpy.commands import validate

config = click.RichHelpConfiguration(
    max_width=80,
    theme = "green2-slim",
    use_markdown=True,
    show_arguments=False,
    style_options_panel_border = "blue",
    style_commands_panel_border = "blue",
    style_option_default= "dim",
    style_deprecated="dim red",
    options_table_column_types = ["opt_long", "opt_short", "help"],
    options_table_help_sections = ["required", "help", "default"]
)

@click.group(options_metavar='', context_settings={"help_option_names" : []} )
@click.rich_config(config)
@click.version_option("0.0.0", prog_name="harpy", hidden = True)
@click.command_panel(
    "Data Processing",
    panel_styles={"border_style": "blue"},
    commands = sorted(["align","demultiplex","qc","snp","sv","impute","phase", "simulate", "assembly", "metassembly"])
)
@click.command_panel(
    "Other Commands",
    panel_styles={"border_style": "magenta"},
    commands = sorted(["convert","deconvolve", "downsample", "template"])
)
@click.command_panel(
    "Troubleshoot",
    panel_styles={"border_style": "yellow"},
    commands = sorted(["view", "resume", "diagnose", "validate", "deps"])
)
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
cli.add_command(environments.deps)
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