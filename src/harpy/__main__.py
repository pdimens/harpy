#!/usr/bin/env python3

import os

try:
    harpypath = '{CONDA_PREFIX}'.format(**os.environ) + "/bin"
except:
    enverror = "\033[1;33mERROR:\033[00m Harpy expects to run from within an active conda environment, but one was not detected."
    print(enverror)
    fix = "\033[1;34mSOLUTION:\033[00m Activate the conda environment Harpy was installed into and run Harpy again."
    print()
    print(fix)
    print(f"\n\033[1mDetails:\033[00m")
    details = "In order to work correctly, Harpy expects several software packages to be available in the PATH, which are provided automatically with Harpy's conda-based installation. It also expects snakefiles, scripts, utilities, etc. to be in the /bin/ folder within that conda environment. If you're seeing this message, no active conda environment was detected upon execution, and Harpy exited as an early failsafe against unexpected runtime errors associated with \"missing\" files and packages."
    print(details)
    exit(1)

from .harpymisc import getnames_err, getnames, vcfcheck
from .extra import extra
from .demultiplex import demultiplex
from .trim import trim
from .align import align
from .variants_snp import snp
from .variants_sv import sv
from .impute import impute
from .phase import phase
import rich_click as click

click.rich_click.USE_MARKDOWN = True
click.rich_click.SHOW_ARGUMENTS = False
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = False
click.rich_click.REQUIRED_SHORT_STRING = ""
click.rich_click.ERRORS_SUGGESTION = "Try the '--help' flag for more information."
click.rich_click.ERRORS_EPILOGUE = "See the documentation: [link=https://pdimens.github.io/harpy/]https://pdimens.github.io/harpy/[/link]"

@click.group(options_metavar='', context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option("0.3", prog_name="Harpy")
def cli():
    """
    ## Harpy haplotagging pipeline
    
    An automated workflow to demultiplex sequences, trim reads, 
    map sequences, call variants, impute genotypes, and phase 
    haplotypes of Haplotagging data. Batteries included.
    
    **demultiplex >> trim >> align >> variants >> impute >> phase**
    
    **Documentation**: [https://pdimens.github.io/harpy/](https://pdimens.github.io/harpy/)
    """
    pass

@click.group(options_metavar='', context_settings=dict(help_option_names=["-h", "--help"]))
def variants():
    """
    Call variants (SNP/SV) from samples

    Provide an additional command `snp` or `sv` to get more information on calling
    those types of variants.
    """
    pass

cli.add_command(extra)
cli.add_command(demultiplex)
cli.add_command(trim)
cli.add_command(align)
cli.add_command(variants)
cli.add_command(impute)
cli.add_command(phase)
variants.add_command(sv)
variants.add_command(snp)

## the modules ##
click.rich_click.OPTION_GROUPS = {
    "harpy demultiplex": [
        {
            "name": "Configuration",
            "options": ["--file", "--samplesheet", "--method"],
        },
        {
            "name": "Other Options",
            "options": ["--threads", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy trim": [
        {
            "name": "Configuration",
            "options": ["--directory", "--max-length", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--threads", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy align": [
        {
            "name": "Configuration",
            "options": ["--genome", "--directory", "--quality-filter", "--method", "--ema-bins", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--threads", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy variants snp": [
        {
            "name": "Configuration",
            "options": ["--genome", "--directory", "--populations", "--ploidy", "--windowsize", "--method", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--threads", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy variants sv": [
        {
            "name": "Configuration",
            "options": ["--genome", "--directory", "--populations", "--method", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--threads", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy impute": [
        {
            "name": "Configuration",
            "options": ["--vcf", "--directory", "--parameters"],
        },
        {
            "name": "Other Options",
            "options": ["--vcf-samples", "--threads", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy phase": [
        {
            "name": "Configuration",
            "options": ["--vcf", "--directory", "--molecule-distance", "--indels", "--prune-threshold", "--ignore-bx", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--vcf-samples", "--threads", "--snakemake", "--quiet", "--help"],
        },
    ]
}

def main():
    cli()

if __name__ == '__main__':
    sys.exit(main())