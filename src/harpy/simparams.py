import os
import sys
from .printfunctions import print_notice
import rich_click as click

genome_conf = """# HARPY GENOME SIMULATION CONFIGURATION FILE #
# ======================================== PLEASE READ ME FIRST ========================================= #
# All of these parameters are optional
# Most simuG parameters are shown here and all with their default values (except random_seed)
# If you aren't interested in simulating a particular type of variant, remove or comment out that section
# If you don't need a particular parameter, remove or comment it
# If you are also using VCF files as inputs, all parameters listed here for those variant types will be ignored
# The bottom of this file has a description of every parameter
# Refer to the Harpy documentation for more details: https://pdimens.github.io/harpy
# ======================================================================================================= #

snp:
    snp_count:
    transition_transversion_ratio: 0.5
indel:
    indel_count:
    insertion_deletion_ratio: 1.0
cnv:
    cnv_count:
    cnv_min_size: 100
    cnv_max_size: 100000
    cnv_max_copy_number: 10
    cnv_gain_loss_ratio:
    tandem_dispersed_ratio: 1
inversion:
    inversion_count:
    inversion_min_size: 1000
    inversion_max_size: 100000
    inversion_breakpoints:
translocation:
    translocation_count:
    translocation_breakpoints:
random_seed: 12343212345
centromere_gff:
gene_gff:
coding_partition_for_snp_simulation:


# ============================================= Parameter descriptions ======================================================= #
## ------------------------------ Randomly simulate Single Nucleotide Polymorphisms (SNP) ----------------------------------- ##
# snp_count: specify the number of SNP variants to be introduced
# transition_transversion_ratio: the transition/transversion ratio used for SNP variants
#     - transition_transversion_ratio: Inf <- only transitions 
#     - transition_transversion_ratio: 0   <- only transversions

## ------------------------------- Randomly simulate  and insertion-deletions (indel) --------------------------------------- ##
# indel_count: the number of indel variants to be introduced
# insertion_deletion_ratio: the insertion/deletion ratio used to simulate indel variants
#    - insertion_deletion_ratio: Inf <- only insertions
#    - insertion_deletion_ratio: 0   <- only deletions

## -------------------------------------- Randomly simulate Copy Number Variants (CNV) ------------------------------------- ##
# cnv_count: number of CNV variants to be introduced
# cnv_min_size: minimum copy number size in base pairs
# cnv_max_size: maximum copy number size in base pairs
# cnv_max_copy_number: maximum number of copies
# cnv_gain_loss_ratio: the relative ratio of DNA gain over DNA loss
#    - cnv_gain_loss_ratio: Inf <- only copy number gain
#    - cnv_gain_loss_ratio: 0   <- only copy number loss
# tandem_dispersed_ratio: expected frequency ratio between tandem and dispersed duplication for CNV variants
#     - tandem_dispersed_ratio: Inf <- only tandem duplications
#     - tandem_dispersed_ratio: 0   <- only dispersed duplications

## --------------------------------------------- Randomly simulate inversions ---------------------------------------------- ##
# inversion_count: number of inversions to be introduced
# inversion_min_size: minimum inversion size in base pairs
# inversion_max_size: maximum inversion size in basepairs
# inversion_breakpoints: GFF3 file of potential breakpoints for triggering inversions

## ------------------------------------------- Randomly simulate translocations -------------------------------------------- ##
# translocation_count: number of translocations to be introduced
# translocation_breakpoints: GFF3 file of the potential breakpoints for triggering translocation

## --------------------------------------------------- Other parameters ---------------------------------------------------- ##
# random_seed: an integer set as the random seed for the simulation
#     - try to set a very large number
# centromere_gff: GFF3 file to specify centromeres for constraining the random CNV, inversion,
#     and translocation simulation. If provided, will prevent those variants types from being
#     simulated in centromeric regions. Potential translocations that will create dicentric
#     rearranged chromosomes will also be disabled. 
# gene_gff: GFF3 file of genes for constraining the random SNP, CNV, inversion, and translocation simulation.
#     - for random CNV, inversion, and translocation simulation, this option will prevent simulated
#       breakpoints falling on the defined genes
#     - for random SNP simulation, this option NEEDS to be used with the option 'coding_partition_for_snp_simulation'
#       to constrain SNP simulations only in noncoding regions, coding regions, 2-fold degenerate (2d) sites or 4-fold degenerate (4d) sites 
# coding_partition_for_snp_simulation: the coding partition used to constrain randomly simulated SNPs 
#     - options: noncoding, coding, 2d, 4d
#     - this option needs to be used together with gene_gff"""

@click.command(no_args_is_help = True, epilog = "read the docs for more information: https://pdimens.github.io/harpy/modules/impute/#parameter-file")
@click.option('-o', '--output', type=str, required = True, help = 'Name of output simulation parameter file')
def simparams(output):
    """
    Create a template parameter file for `harpy simulate genome`

    With this command you can create a template parameter
    file necessary for simulating random variants in a genome via `harpy simulate`
    The resulting file will have default `simuG` values and should be modified appropriately.
    """
    if os.path.exists(output):
        overwrite = input(f"File {output} already exists, overwrite (no|yes)?  ").lower()
        if overwrite not in ["yes", "y"]:
            click.echo("Please suggest a different name for the output file")
            exit(0)

    with open(output, "w+") as f:
        f.writelines(genome_conf)

    print_notice(
        f"Created simulation parameter file: {output}\n" +
        "Modify the file as needed, details are provided inside the file."
    )

