#!/usr/bin/env python3

import rich_click as click
from . import align
from . import diagnose, resume, view
from . import deconvolve
from . import demultiplex
from . import container
from . import hpc
from . import impute
from . import assembly, metassembly
from . import qc
from . import phase
from . import preflight
from . import simulate_linkedreads, simulate_variants
from . import snp
from . import sv
from .popgroup import popgroup
from . import downsample
from .imputeparams import imputeparams

click.rich_click.USE_MARKDOWN = True
click.rich_click.SHOW_ARGUMENTS = False
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = False
click.rich_click.MAX_WIDTH = 75
click.rich_click.REQUIRED_SHORT_STRING = ""
click.rich_click.ERRORS_SUGGESTION = "Try the '--help' flag for more information."
click.rich_click.ERRORS_EPILOGUE = "Documentation: [link=https://pdimens.github.io/harpy/]https://pdimens.github.io/harpy/[/link]"

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

## unify simulate commands
@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def simulate():
    """
    Simulate variants or linked-reads from a genome

    To simulate genomic variants, provide an additional subcommand {`snpindel`,`inversion`,`cnv`,`translocation`} 
    to get more information about that workflow. The variant simulator (`simuG`) can only simulate
    one type of variant at a time, so you may need to run it a few times if you want multiple variant types.
    Use `simulate linkedreads` to simulate haplotag linked-reads from a diploid genome, which you can create by simulating
    genomic variants.
    """

simulate_commandstring = {
    "harpy simulate": [
        {
            "name": "Linked Read Sequences",
            "commands": ["linkedreads"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Genomic Variants",
            "commands": ["cnv", "inversion", "snpindel", "translocation"],
            "panel_styles": {"border_style": "green"}
        }
    ]
}

simulate.add_command(simulate_linkedreads.linkedreads)
simulate.add_command(simulate_variants.snpindel)
simulate.add_command(simulate_variants.inversion)
simulate.add_command(simulate_variants.cnv)
simulate.add_command(simulate_variants.translocation)

# main program
cli.add_command(downsample.downsample)
cli.add_command(popgroup)
cli.add_command(imputeparams)
cli.add_command(view.view)
cli.add_command(preflight.preflight)
cli.add_command(demultiplex.demultiplex)
cli.add_command(qc.qc)
cli.add_command(align.align)
cli.add_command(snp.snp)
cli.add_command(sv.sv)
cli.add_command(impute.impute)
cli.add_command(phase.phase)
cli.add_command(simulate)
cli.add_command(container.containerize)
cli.add_command(hpc.hpc)
cli.add_command(resume.resume)
cli.add_command(deconvolve.deconvolve)
cli.add_command(metassembly.metassembly)
cli.add_command(assembly.assembly)
cli.add_command(diagnose.diagnose)
click.rich_click.COMMAND_GROUPS = {
    "harpy":
        [
            {
                "name": "Workflows",
                "commands": sorted(["demultiplex","qc", "align","snp","sv","impute","phase", "simulate", "assembly", "metassembly"]),
                "panel_styles": {"border_style": "blue"}
            },
            {
                "name": "Other Commands",
                "commands": sorted(["deconvolve", "downsample", "hpc", "imputeparams", "popgroup"]),
                "panel_styles": {"border_style": "dim green"}
            },
            {
                "name": "Troubleshoot",
                "commands": sorted(["view", "resume", "diagnose", "preflight"]),
                "panel_styles": {"border_style": "dim"}
            }
        ],
 } | simulate_commandstring | hpc.docstring

click.rich_click.OPTIONS_PANEL_TITLE = None
for i in [align, deconvolve, downsample, demultiplex, impute, phase, preflight, qc, simulate_linkedreads, simulate_variants, snp, sv, assembly, metassembly]:
    click.rich_click.OPTION_GROUPS |= i.docstring
