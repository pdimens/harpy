#!/usr/bin/env python3

import rich_click as click

from harpy import __version__
from harpy.common.lazy_group import LazyGroup, LazySpec

class RichLazyGroup(LazyGroup, click.RichGroup):
    """Lazily-importing group with rich-click help formatting"""

# Subcommands are imported only when run. The help text (first paragraph of each command's docstring)
# and hidden flag are repeated here so the top-level `--help` and shell completion can list them
# without importing anything. test/test_lazy_cli.py checks that this stays in sync with the commands.
COMMANDS = {
    "align":        LazySpec("harpy.commands.align:align", "Align sequences to a reference genome"),
    "assembly":     LazySpec("harpy.commands.assembly:assembly", "Assemble linked reads into a genome"),
    "completion":   LazySpec("harpy.commands.completion:completion", "Print the shell completion script", hidden = True),
    "containerize": LazySpec("harpy.commands.environments:containerize", "Configure the harpy container", hidden = True),
    "deconvolve":   LazySpec("harpy.commands.deconvolve:deconvolve", "Resolve barcode sharing in unrelated molecules"),
    "deps":         LazySpec("harpy.commands.environments:deps", "Locally install workflow dependencies"),
    "diagnose":     LazySpec("harpy.commands.diagnose:diagnose", "Attempt to resolve workflow errors"),
    "impute":       LazySpec("harpy.commands.impute:impute", "Impute variant genotypes from alignments"),
    "metassembly":  LazySpec("harpy.commands.assembly:metassembly", "Assemble linked reads into a metagenome"),
    "phase":        LazySpec("harpy.commands.phase:phase", "Phase SNPs or alignments"),
    "preprocess":   LazySpec("harpy.commands.preprocess:preprocess", "Remove inline barcodes from raw FASTQs"),
    "qc":           LazySpec("harpy.commands.qc:qc", "FASTQ adapter removal, quality filtering, etc."),
    "report":       LazySpec("harpy.commands.report:report", "Render ipynb reports"),
    "resume":       LazySpec("harpy.commands.resume:resume", "Continue an incomplete Harpy workflow"),
    "snp":          LazySpec("harpy.commands.snp:snp", "Call SNPs and small indels from alignments"),
    "sv":           LazySpec("harpy.commands.sv:sv", "Call inversions, deletions, and duplications from alignments"),
    "template":     LazySpec("harpy.commands.template:template", "Create files and HPC configs for workflows"),
    "validate":     LazySpec("harpy.commands.validate:validate", "File format checks for linked-read data"),
    "view":         LazySpec("harpy.commands.view:view", "View a workflow's components"),
}

config = click.RichHelpConfiguration(
    max_width=80,
    theme = "green2-slim",
    use_markdown=True,
    show_arguments=False,
    style_options_panel_border = "blue",
    style_commands_panel_border = "blue",
    style_option_default= "dim",
    style_deprecated="dim red",
    style_errors_panel_border = "yellow",
    errors_panel_title = "Usage Error",
    options_table_column_types = ["opt_long", "opt_short", "help"],
    options_table_help_sections = ["required", "help", "default"]
)

@click.group(cls = RichLazyGroup, lazy_commands = COMMANDS, options_metavar='')
@click.rich_config(config)
@click.version_option(__version__, prog_name="harpy", hidden = True)
@click.command_panel(
    "Data Processing",
    panel_styles={"border_style": "blue"},
    commands = sorted(["align", "deconvolve", "preprocess","qc","snp","sv","impute","phase", "assembly", "metassembly"])
)
@click.command_panel(
    "Other Commands",
    panel_styles={"border_style": "magenta"},
    commands = sorted(["report", "template"])
)
@click.command_panel(
    "Troubleshoot",
    panel_styles={"border_style": "yellow"},
    commands = sorted(["view", "resume", "diagnose", "validate", "deps"])
)
@click.help_option('--help', hidden = True)
def cli():
    """
    Automated workflows for linked-read data
    to go from raw data to genotypes (or phased haplotypes).
    Batteries included.

    **preprocess >> qc >> align >> snp >> impute >> phase >> sv**

    **Documentation**: [https://pdimens.github.io/harpy/](https://pdimens.github.io/harpy/)
    """

if __name__ == "__main__":
    cli()
