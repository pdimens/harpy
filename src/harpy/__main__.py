#!/usr/bin/env python3

import os
import sys

try:
    harpypath = '{CONDA_PREFIX}'.format(**os.environ) + "/bin"
except:
    enverror = "\033[1;33mERROR:\033[00m Harpy expects to run from within an active conda environment, but one was not detected."
    print(enverror, file = sys.stderr)
    fix = "\033[1;34mSOLUTION:\033[00m Activate the conda environment Harpy was installed into and run Harpy again."
    print()
    print(fix, file = sys.stderr)
    print(f"\n\033[1mDetails:\033[00m")
    details = "In order to work correctly, Harpy expects several software packages to be available in the PATH, which are provided automatically with Harpy's conda-based installation. It also expects snakefiles, scripts, utilities, etc. to be in the /bin/ folder within that conda environment. If you're seeing this message, no active conda environment was detected upon execution, and Harpy exited as an early failsafe against unexpected runtime errors associated with \"missing\" files and packages."
    print(details, file = sys.stderr)
    exit(1)

from . import align
from . import demultiplex
from . import impute
from . import phase
from . import preflight
from . import qc
from . import simulate
from . import snp
from . import sv
from .popgroups import popgroup
from .stitchparams import stitchparams
from .conda_deps import generate_conda_deps
import rich_click as click
import subprocess

click.rich_click.USE_MARKDOWN = True
click.rich_click.SHOW_ARGUMENTS = False
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = False
click.rich_click.MAX_WIDTH = 75
click.rich_click.REQUIRED_SHORT_STRING = ""
click.rich_click.ERRORS_SUGGESTION = "Try the '--help' flag for more information."
click.rich_click.ERRORS_EPILOGUE = "See the documentation: [link=https://pdimens.github.io/harpy/]https://pdimens.github.io/harpy/[/link]"

@click.group(options_metavar='', context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option("1.0.0", prog_name="Harpy")
def cli():
    """
    ## Harpy haplotagging pipeline
    
    An automated workflow to demultiplex sequences, trim and qc reads, 
    map sequences, call variants, impute genotypes, and phase 
    haplotypes of Haplotagging data. Batteries included.
    
    **demultiplex >> qc >> align >> snp >> impute >> phase >> sv**
    
    **Documentation**: [https://pdimens.github.io/harpy/](https://pdimens.github.io/harpy/)
    """
    pass

# main program
#cli.add_command(hpc)
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
# demultiplex submodules
demultiplex.demultiplex.add_command(demultiplex.gen1)
# preflight submodules
preflight.preflight.add_command(preflight.fastq)
preflight.preflight.add_command(preflight.bam)
# align submodules
align.align.add_command(align.bwa)
align.align.add_command(align.ema)
align.align.add_command(align.minimap)
# snp submodules
snp.snp.add_command(snp.mpileup)
snp.snp.add_command(snp.freebayes)
# sv submodules
sv.sv.add_command(sv.leviathan)
sv.sv.add_command(sv.naibr)
# simulate submodules
simulate.simulate.add_command(simulate.snpindel)
simulate.simulate.add_command(simulate.inversion)
simulate.simulate.add_command(simulate.cnv)
simulate.simulate.add_command(simulate.translocation)
simulate.simulate.add_command(simulate.linkedreads)
#simulate.add_command(reads)

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
                "commands": ["preflight", "popgroup", "stitchparams"]
            }
        ],
    "harpy simulate":
        [
            {
                "name": "Linked Read Sequences",
                "commands": ["linkedreads"],
            },
            {
                "name": "Genomic Variants",
                "commands": ["snpindel","inversion", "cnv", "translocation"],
            }
        ]
}

click.rich_click.OPTION_GROUPS = demultiplex.docstring | preflight.docstring | qc.docstring | align.docstring | snp.docstring | sv.docstring | impute.docstring | phase.docstring | simulate.docstring

def main():
    try:
        workflow = cli(standalone_mode = False)
        if workflow == 0:
            return 0
        elif workflow is not None:
            generate_conda_deps()
            _module = subprocess.run(workflow, check = True)
            return _module.returncode
    except:
        sys.exit(1)

if __name__ == '__main__':
    sys.exit(main())