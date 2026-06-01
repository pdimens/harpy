import sys

import click

@click.command(no_args_is_help = False)
@click.option('-m', '--me-seq', type=str, default = "AGATGTGTATAAGAGACAG", show_default = True, help = "The ME sequence used in tagmentation")
@click.help_option('--help', hidden = True)
def known_adapters(me_seq):
    """
    INTERNAL USE- Writes a fasta file of common illumina adapters that might appear in the GIH Nextera prep, the one used by Cornell GIH for tagmentation, and the ME sequence
    for use in QC adapter trimming. Writes to stdout. Most of these were derived from fastp (https://github.com/OpenGene/fastp)
    """
    fa = f""">Nextera GIH
CTGTCTCTTATACACATCT
>ME Sequence
{me_seq}
>I7_Nextera_Transposase_1 | >Trans2_rc | >I7_Nextera_Transposase_1 | >Trans2_rc
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
>I7_Nextera_Transposase_2
CTGTCTCTTATACACATCTCTGAGCGGGCTGGCAAGGC
>I5_Nextera_Transposase_2
CTGTCTCTTATACACATCTCTGATGGCGCGAGGGAGGC
>I5_Nextera_Transposase_1 | >Trans1_rc | >I5_Nextera_Transposase_1 | >Trans1_rc
CTGTCTCTTATACACATCTGACGCTGCCGACGA
"""
    sys.stdout.write(fa)