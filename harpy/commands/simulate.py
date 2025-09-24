"""Unify the harpy simulate commands"""

import rich_click as click
from . import simulate_variants

@click.group(options_metavar='', context_settings={"help_option_names" : []})
@click.command_panel("Linked Read Sequences", commands = ["linkedreads"])
@click.command_panel("Genomic Variants", commands = ["cnv", "inversion", "snpindel", "translocation"])
def simulate():
    """
    Simulate genomic variants

    To simulate genomic variants, provide an additional subcommand {`snpindel`,`inversion`,`cnv`,`translocation`} 
    to get more information about that workflow. The variant simulator (`simuG`) can only simulate
    one type of variant at a time, so you may need to run it a few times if you want multiple variant types.
    """

import sys
import rich_click as click

@click.command(deprecated=True, no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/mimick/")
def linkedreads():
    """
    Create linked reads using genome haplotypes

    This workflow was a very thin veneer to use `Mimick` with mutations (SNPs and indels) deactivated. So, 
    we now recommend using Mimick directly to achieve this by installing `Mimick` from Bioconda and using it directly.
    """
    sys.exit(1)

simulate.add_command(linkedreads)
simulate.add_command(simulate_variants.snpindel)
simulate.add_command(simulate_variants.inversion)
simulate.add_command(simulate_variants.cnv)
simulate.add_command(simulate_variants.translocation)