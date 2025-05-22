"""Downsample a fastq/bam file by barcodes"""

import os
import re
import sys
from turtle import down
import rich_click as click
from ._cli_types_generic import convert_to_int, SnakemakeParams, HPCProfile
from ._launch import launch_snakemake
from ._misc import fetch_rule, instantiate_dir, setup_snakemake, write_workflow_config
from ._printing import workflow_info

docstring = {
    "harpy downsample": [
        {
            "name": "Parameters",
            "options": sorted(["--barcode-tag", "--downsample", "--invalid", "--random-seed", "--prefix"]),
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--hpc", "--output-dir", "--quiet", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
}

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/downsample")
@click.option('-b', '--barcode-tag', type = str, default = "BX", show_default = True, help = 'Tag that contains the barcode')
@click.option('-d', '--downsample', type = click.FloatRange(min = 0.0000001), help = 'Number/fraction of barcodes to retain')
@click.option('-i', '--invalid', default = 1, show_default = True, type=click.FloatRange(min=0,max=1), help = "Proportion of invalid barcodes to sample")
@click.option('-p', '--prefix', type = click.Path(exists = False), default = "downsampled", show_default = True, help = 'Prefix for output file(s)')
@click.option('--random-seed', type = click.IntRange(min = 1), help = "Random seed for sampling")
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Downsample", show_default=True,  help = 'Output directory name')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(1,999, clamp = True), help = 'Number of threads to use')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = "0", type = click.Choice(["0", "1", "2"]), callback = convert_to_int, help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--setup-only',  is_flag = True, hidden = True, show_default = True, default = False, help = 'Setup the workflow and exit')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('input', required=True, type=click.Path(exists=True, readable=True, dir_okay=False, resolve_path=True), nargs=-1)
def downsample(input, invalid, output_dir, prefix, barcode_tag, downsample, random_seed, hpc, setup_only, snakemake, threads, quiet):
    """
    Downsample data by barcode
    
    Downsamples FASTQ or BAM file(s) by barcode in the `BX` (default) tag to keep all reads
    containing `-d` randomly sampled barcodes. Downsamples by keeping `-d` random barcodes if `d >= 1`,
    otherwise, samples a fraction of the total barcodes if `0 < d < 1` (e.g. `-d 0.5` retains 50% of all barcodes).
     The SAM tag TYPE is expected to be `Z`, i.e. `--barcode-type BX` will search through the `BX:Z:` tags. 
    Use `--invalid/-i` to specify the proportion of invalid barcodes to consider for sampling. Input can be:
    - one BAM file
    - two FASTQ files (R1 and R2 of a paired-end read set)
    
    | `--invalid` | effect                                           |
    |:---:|:---------------------------------------------------|
    | `0` | removes all invalid barcodes from the sampling pool |
    | `1` | adds all invalid barcodes to the sampling pool |
    | 0<`i`<1| keeps `i` proprotion of invalids in the sampling pool |
    """
    workflow = "downsample"
    workflowdir,sm_log = instantiate_dir(output_dir, workflow)
    ## checks and validations ##
    if len(input) > 2:
        raise click.BadParameter('inputs must be 1 BAM file or 2 FASTQ files.', param_hint="INPUT")
    if len(input) == 1:
        if not input[0].lower().endswith(".bam"):
            raise click.BadParameter('inputs must be 1 BAM (.bam) file or 2 FASTQ (.fastq|.fq) files. The FASTQ files can be gzipped.', param_hint="INPUT")
    if len(input) == 2:
        if input[0] == input[1]:
            raise click.BadParameter('the two input files cannot be identical', param_hint="INPUT")
        re_ext = re.compile(r"\.(fq|fastq)(?:\.gz)?$", re.IGNORECASE)
        for i in input:
            if not re.search(re_ext, i):
                raise click.BadParameter('inputs must be 1 BAM (.bam) file or 2 FASTQ (.fastq|.fq) files. The FASTQ files can be gzipped.', param_hint="INPUT")            
    if len(barcode_tag) != 2:
        raise click.BadParameter('The barcode tag must be 2 chracters from the English alphabet (A-Z)', param_hint="--barcode-tag/-b")            


    ## setup workflow ##
    command,command_rel = setup_snakemake(
        workflow,
        "conda",
        output_dir,
        threads,
        hpc if hpc else None,
        snakemake if snakemake else None
    )

    fetch_rule(workflowdir, "downsample.smk")

    configs = {
        "workflow": workflow,
        "prefix" :  prefix,
        "downsample" :  int(downsample) if downsample >= 1 else downsample,
        "barcode-tag" : barcode_tag.upper(),
        "invalid_proportion" : invalid,       
        **({"random_seed" : random_seed} if random_seed else {}),
        "snakemake" : {
            "log" : sm_log,
            "absolute": command,
            "relative": command_rel
        },
        "inputs": input
    }

    write_workflow_config(configs, output_dir)
    if setup_only:
        sys.exit(0)

    start_text = workflow_info(
        ("Output Folder:", os.path.basename(output_dir) + "/"),
        ("Downsample barcodes:" if downsample >= 1 else "Downsample fraction:", f"{int(downsample) if downsample >= 1 else downsample}"),
        ("Invalid Proportion:", invalid),
        ("Workflow Log:", sm_log.replace(f"{output_dir}/", "") + "[dim].gz")
    )
    launch_snakemake(command_rel, workflow, start_text, output_dir, sm_log, quiet, f"workflow/downsample.summary")
