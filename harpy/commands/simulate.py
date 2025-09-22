"""Unify the harpy simulate commands"""

import rich_click as click
from . import simulate_linkedreads
from . import simulate_variants

@click.group(options_metavar='', context_settings={"help_option_names" : []})
@click.command_panel("Linked Read Sequences", commands = ["linkedreads"])
@click.command_panel("Genomic Variants", commands = ["cnv", "inversion", "snpindel", "translocation"])
def simulate():
    """
    Simulate genomic variants or linked reads

    To simulate genomic variants, provide an additional subcommand {`snpindel`,`inversion`,`cnv`,`translocation`} 
    to get more information about that workflow. The variant simulator (`simuG`) can only simulate
    one type of variant at a time, so you may need to run it a few times if you want multiple variant types.
    Use `simulate linkedreads` to simulate haplotagging linked reads from a diploid genome, which you can create by simulating
    genomic variants.
    """

simulate.add_command(simulate_linkedreads.linkedreads)
simulate.add_command(simulate_variants.snpindel)
simulate.add_command(simulate_variants.inversion)
simulate.add_command(simulate_variants.cnv)
simulate.add_command(simulate_variants.translocation)