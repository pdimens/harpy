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
from .stitchparams import stitchparams
from .hpc import hpc
from .demultiplex import gen1
from .preflight import bam, fastq
from .qc import qc
from .align import bwa, ema
from .snp import freebayes, mpileup
from .sv import leviathan, naibr
from .impute import impute
from .phase import phase
import rich_click as click

click.rich_click.USE_MARKDOWN = True
click.rich_click.SHOW_ARGUMENTS = False
click.rich_click.SHOW_METAVARS_COLUMN = False
click.rich_click.APPEND_METAVARS_HELP = False
click.rich_click.MAX_WIDTH = 75
click.rich_click.REQUIRED_SHORT_STRING = ""
click.rich_click.ERRORS_SUGGESTION = "Try the '--help' flag for more information."
click.rich_click.ERRORS_EPILOGUE = "See the documentation: [link=https://pdimens.github.io/harpy/]https://pdimens.github.io/harpy/[/link]"

@click.group(options_metavar='', context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option("0.9.0", prog_name="Harpy")
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

    **Aligners**
    - `bwa`: uses BWA MEM to align reads, retaining BX tags in the alignments
    - `ema`: uses the BX barcode-aware EMA aligner

    Provide an additional subcommand `bwa` or `ema` to get more information on using
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
    so you may need to run this module again on the resulting genome.
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
                "name": "Genomic Variants",
                "commands": ["snpindel","inversion", "cnv", "translocation"],
            },
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
            "name": "Module Parameters",
            "options": ["--samplesheet"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy qc": [
        {
            "name": "Module Parameters",
            "options": ["--max-length", "--ignore-adapters", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy align bwa": [
        {
            "name": "Module Parameters",
            "options": ["--genome", "--quality-filter", "--molecule-distance", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy align ema": [
        {
            "name": "Module Parameters",
            "options": ["--platform", "--whitelist", "--genome", "--quality-filter", "--ema-bins", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy snp mpileup": [
        {
            "name": "Module Parameters",
            "options": ["--genome", "--populations", "--ploidy", "--windowsize", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy snp freebayes": [
        {
            "name": "Module Parameters",
            "options": ["--genome", "--populations", "--ploidy", "--windowsize", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy sv leviathan": [
        {
            "name": "Module Parameters",
            "options": ["--genome", "--populations", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy sv naibr": [
        {
            "name": "Module Parameters",
            "options": ["--genome", "--vcf", "--molecule-distance", "--populations", "--extra-params"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy impute": [
        {
            "name": "Module Parameters",
            "options": ["--vcf", "--parameters", "--extra-params", "--vcf-samples"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--snakemake", "--quiet", "--help"],
        },
    ],
    "harpy phase": [
        {
            "name": "Module Parameters",
            "options": ["--vcf", "--molecule-distance", "--genome", "--prune-threshold", "--ignore-bx", "--extra-params", "--vcf-samples"],
        },
        {
            "name": "Other Options",
            "options": ["--output-dir", "--threads", "--skipreports", "--snakemake", "--quiet", "--help"],
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
            "options": ["--output-dir", "--heterozygosity", "--randomseed", "--snakemake", "--quiet", "--help"],
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
            "options": ["--output-dir", "--heterozygosity", "--randomseed", "--snakemake", "--quiet", "--help"],
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
            "options": ["--output-dir", "--heterozygosity", "--randomseed", "--snakemake", "--quiet", "--help"],
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
            "options": ["--output-dir","--heterozygosity", "--randomseed", "--snakemake", "--quiet", "--help"],
        },
    ],
}


def main():
    cli()

if __name__ == '__main__':
    sys.exit(main())