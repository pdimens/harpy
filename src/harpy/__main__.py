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

from .popgroups import popgroup
#from .simulatelinkedreads import reads
from .simulatevariants import snpindel, inversion, cnv, translocation
from .simulatelinkedreads import linkedreads
from .stitchparams import stitchparams
from .hpc import hpc
from .demultiplex import gen1
from .preflight import bam, fastq
from .qc import qc
from .align import bwa, ema, minimap
from .snp import freebayes, mpileup
from .sv import leviathan, naibr
from .impute import impute
from .phase import phase
import rich_click as click
import subprocess
from .helperfunctions import generate_conda_deps

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

@click.group(options_metavar='', context_settings=dict(help_option_names=["-h", "--help"]))
def demultiplex():
    """
    Demultiplex haplotagged FASTQ files

    Check that you are using the correct haplotag method/technology, since the different
    barcoding approaches have very different demultiplexing strategies.

    **Haplotag Technologies**
    - `gen1`: the original haplotagging barcode strategy developed by Meier _et al._ (2021)
    """

@click.group(options_metavar='', context_settings=dict(help_option_names=["-h", "--help"]))
def align():
    """
    Align sample sequences to a reference genome

    The three available aligners all retain the linked-read barcode information in the
    resulting output, however `EMA` is the only aligner to use the barcode information
    to facilitate the aligning process and can be prohibitively slow. The `minimap2`
    aligner is the fastest of the three and is comparable in accuracy to `bwa` for
    sequences >100bp.

    **Aligners**
    - `bwa`: uses BWA MEM to align reads (fast)
    - `ema`: uses the BX barcode-aware EMA aligner (very slow)
    - `minimap`: uses minimap2 to align reads (ultra fast)

    Provide an additional subcommand `bwa`, `ema`, or `minimap` to get more information on using
    those aligners.
    """
    pass

@click.group(options_metavar='', context_settings=dict(help_option_names=["-h", "--help"]))
def snp():
    """
    Call SNPs and small indels
    
    **Variant Callers**
    - `mpileup`: call variants using bcftools mpileup
    - `freebayes`: call variants using freebayes

    Provide an additional subcommand `mpileup` or `freebayes` to get more information on using
    those variant callers. They are both robust variant callers and neither is recommended over the other.
    """
    pass

@click.group(options_metavar='', context_settings=dict(help_option_names=["-h", "--help"]))
def sv():
    """
    Call large structural variants
 
    **Structural Variant Callers**
    - `naibr`: calls inversions, duplicates, deletions
    - `leviathan`: calls inversions, duplicates, deletions, misc breakends

    Provide an additional subcommand `leviathan` or `naibr` to get more information on using
    those variant callers. NAIBR tends to call variants better, but requires more user preprocessing.
    """
    pass

@click.group(options_metavar='', context_settings=dict(help_option_names=["-h", "--help"]))
def preflight():
    """
    Run file format checks on haplotag data

    This is useful to make sure your input files are formatted correctly for the processing pipeline 
    before you are surprised by errors hours into an analysis. Provide an additional command `fastq`
    or `bam` to see more information and options.
    """
    pass

@click.group(options_metavar='', context_settings=dict(help_option_names=["-h", "--help"]))
def simulate():
    """
    Simulate variants or linked reads from a genome

    To simulate genomic variants, provide an additional subcommand {`snpindel`,`inversion`,`cnv`,`translocation`} 
    to get more information about that workflow. The limitations of the simulator
    (`simuG`) are such that you may simulate only one type of variant at a time,
    so you may need to run this module again on the resulting genome. Use `simulate linkedreads`
    to simulate haplotag linked-reads from a diploid genome, which you can create by simulating
    genomic variants.
    """
    pass

# main program
#cli.add_command(hpc)
cli.add_command(popgroup)
cli.add_command(stitchparams)
cli.add_command(preflight)
cli.add_command(demultiplex)
cli.add_command(qc)
cli.add_command(align)
cli.add_command(snp)
cli.add_command(sv)
cli.add_command(impute)
cli.add_command(phase)
cli.add_command(simulate)
# demultiplex submodules
demultiplex.add_command(gen1)
# preflight submodules
preflight.add_command(fastq)
preflight.add_command(bam)
# align submodules
align.add_command(bwa)
align.add_command(ema)
align.add_command(minimap)
# snp submodules
snp.add_command(mpileup)
snp.add_command(freebayes)
# sv submodules
sv.add_command(leviathan)
sv.add_command(naibr)
# simulate submodules
simulate.add_command(snpindel)
simulate.add_command(inversion)
simulate.add_command(cnv)
simulate.add_command(translocation)
simulate.add_command(linkedreads)
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

click.rich_click.OPTION_GROUPS = {
    "harpy preflight bam": [
        {
            "name": "Options",
            "options": ["--output-dir", "--threads", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy preflight fastq": [
        {
            "name": "Options",
            "options": ["--output-dir", "--threads", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy demultiplex gen1": [
        {
            "name": "Parameters",
            "options": ["--schema"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy qc": [
        {
            "name": "Parameters",
            "options": ["--max-length", "--ignore-adapters", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy align bwa": [
        {
            "name": "Parameters",
            "options": ["--genome", "--quality-filter", "--molecule-distance", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy align ema": [
        {
            "name": "Parameters",
            "options": ["--platform", "--whitelist", "--genome", "--quality-filter", "--ema-bins", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy align minimap": [
        {
            "name": "Parameters",
            "options": ["--genome", "--quality-filter", "--molecule-distance", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy snp mpileup": [
        {
            "name": "Parameters",
            "options": ["--genome", "--populations", "--ploidy", "--regions", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy snp freebayes": [
        {
            "name": "Parameters",
            "options": ["--genome", "--populations", "--ploidy", "--regions", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy sv leviathan": [
        {
            "name": "Parameters",
            "options": ["--genome", "--min-sv", "--min-barcodes", "--populations", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy sv naibr": [
        {
            "name": "Module Parameters",
            "options": ["--genome", "--vcf", "--min-sv", "--min-barcodes", "--molecule-distance", "--populations", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy impute": [
        {
            "name": "Parameters",
            "options": ["--vcf", "--parameters", "--extra-params", "--vcf-samples"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy phase": [
        {
            "name": "Parameters",
            "options": ["--vcf", "--molecule-distance", "--genome", "--prune-threshold", "--ignore-bx", "--extra-params", "--vcf-samples"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--snakemake", "--quiet", "--help"],
        },     
    ],
    "harpy simulate linkedreads": [
        {
            "name": "Parameters",
            "options": ["--barcodes", "--read-pairs", "--outer-distance", "--distance-sd", "--molecule-length", "--partitions", "--molecules-per"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--snakemake", "--quiet", "--help"],
        },     
    ],
    "harpy simulate snpindel": [
        {
            "name": "Known Variants",
            "options": ["--snp-vcf", "--indel-vcf"],
        },
        {
            "name": "Random Variants",
            "options": ["--snp-count", "--indel-count", "--titv-ratio", "--indel-ratio", "--snp-gene-constraints", "--genes", "--centromeres", "--exclude-chr"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--prefix", "--heterozygosity", "--randomseed", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy simulate inversion": [
        {
            "name": "Known Variants",
            "options": ["--vcf"],
        },
        {
            "name": "Random Variants",
            "options": ["--count", "--min-size", "--max-size", "--genes", "--centromeres", "--exclude-chr"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--prefix", "--heterozygosity", "--randomseed", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy simulate cnv": [
        {
            "name": "Known Variants",
            "options": ["--vcf"],
        },
        {
            "name": "Random Variants",
            "options": ["--count", "--min-size", "--max-size", "--max-copy", "--dup-ratio", "--gain-ratio", "--genes", "--centromeres", "--exclude-chr"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--prefix", "--heterozygosity", "--randomseed", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy simulate translocation": [
        {
            "name": "Known Variants",
            "options": ["--vcf"],
        },
        {
            "name": "Random Variants",
            "options": ["--count", "--genes", "--centromeres", "--exclude-chr"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir","--prefix","--heterozygosity", "--randomseed", "--snakemake", "--quiet", "--help"],
        },
    ],
}


def main():
    try:
        workflow = cli(standalone_mode = False)
        if workflow == 0:
            return 0
        elif workflow is not None:
            generate_conda_deps()
            _module = subprocess.run(workflow)
            return _module.returncode
    except:
        sys.exit(1)

if __name__ == '__main__':
    sys.exit(main())