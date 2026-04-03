import os
import sys
import click

@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/preprocess/")
@click.argument('input', required = True, type=click.File())
@click.help_option('--help', hidden = True)
def parse_phaseblocks(input):
    """
    Summarize a HapCut2 phase block file
    
    Writes to stdout.
    """
    samplename = os.path.basename(input.name).replace(".blocks", "")

    FIRST_LINE = True
    for line in input:
        lsplit = line.split()
        if lsplit[0] == "BLOCK:":
            n_snp     = int(lsplit[6])
            len_block = int(lsplit[8])
            FIRST_LINE = True
        else:
            if FIRST_LINE:
                pos_start = int(lsplit[4])
                contig    = lsplit[3]
                FIRST_LINE = False
                sys.stdout.write(f"{samplename}\t{contig}\t{n_snp}\t{pos_start}\t{len_block}\n")
