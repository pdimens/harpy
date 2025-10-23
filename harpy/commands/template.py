"""Harpy module to create a sample grouping file"""

import os
import re
import sys
import glob
import rich_click as click
from harpy.common.printing import print_error, print_notice, CONSOLE
from harpy.common.file_ops import fetch_template

@click.group(context_settings={"help_option_names" : []})
@click.command_panel("Input Files", panel_styles={"border_style": "blue"})
@click.command_panel("HPC Configurations", panel_styles={"border_style": "blue"})
def template():
    """
    Create files and HPC configs for workflows

    All commands write to `stdout`. Use hpc-* and impute without arguments.
    """

@click.command(panel = "Input Files", no_args_is_help=True, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/snp/#sample-grouping-file")
@click.argument('inputdir', required=True, type=click.Path(exists=True, file_okay=False))
def groupings(inputdir):
    """
    Create a template sample-grouping file

    This command generates a sample grouping file, like the kind optional for variant calling.
    Provide the input fastq/bam directory at the end of the command. Writes to `stdout`.
    **Note** that Harpy cannot reliably infer populations from filenames, therefore all samples
    will be assigned to `pop1`. Please modify this file with appropriate population
    designations.
    """
    try:
        samplenames = set()
        re_ext = re.compile(r"\.(bam|sam)$", re.IGNORECASE)
        for i in os.listdir(inputdir):
            if i.lower().endswith(".bam") or i.lower().endswith(".sam"):
                samplenames.add(re_ext.sub("", os.path.basename(i)))
        if len(samplenames) < 1:
            raise ValueError
    except ValueError:
       full_flist = [i for i in glob.iglob(f"{inputdir}/*") if not os.path.isdir(i)]
       r = re.compile(r".*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
       full_fqlist = list(filter(r.match, full_flist))
       fqlist = [os.path.basename(i) for i in full_fqlist]
       bn_r = r"[\.\_][RF](?:[12])?(?:\_00[1-9])*\.f(?:ast)?q(?:\.gz)?$"
       if len(fqlist) == 0:
            print_error(
                "no files found",
                f"No [bold]FASTQ[/] or [bold]BAM[/] files were detected in [blue]{inputdir}[/]",
                "Check that [bold]FASTQ[/] file endings conform to [green].[/][[green]F[/][dim]|[/dim][green]R1[/]][green].[/][[green]fastq[/][dim]|[/dim][green]fq[/]][green].gz[/]\nCheck that [bold]BAM[/] files end in [green].bam[/]\nRead the documentation for details: https://pdimens.github.io/harpy/haplotagdata/#naming-conventions"
           )
       samplenames = set([re.sub(bn_r, "", i, flags = re.IGNORECASE) for i in fqlist])

    CONSOLE.print(f"\n[bold]{len(samplenames)}[/] samples detected in [blue]{inputdir}[/]\n")
    for i in samplenames:
        _ = sys.stdout.write(f'{i}\tpop1\n')
    print_notice("Please review the resulting file, as all samples have been grouped into a single population")

import rich_click as click
from harpy.common.file_ops import fetch_template
from harpy.common.system_ops import package_absent

@click.command(panel = "HPC Configurations")
def hpc_generic():
    """
    Create a template config for a generic scheduler
    
    This command creates a configuration for a generic HPC scheduler. Writes to `stdout`.
    You will also need to install `snakemake-executor-plugin-cluster-generic` for the HPC job submission to work.
    """
    fetch_template("hpc-generic.yaml")
    package_absent("snakemake-executor-plugin-cluster-generic")

#def hpc_generic2():
#    pass
    #executor: cluster-generic
    #cluster-generic-submit-cmd:
    #  mkdir -p results/slurm_logs/{rule} &&
    #  sbatch
    #    --partition=regular
    #    --account=username
    #    --cpus-per-task={threads}
    #    --mem-per-cpu={resources.mem_per_cpu}
    #    --time={resources.time}
    #    --job-name=harpy-{rule}-{wildcards}
    #    --output=results/slurm_logs/{rule}/{rule}-%j-{wildcards}.out
    #    --error=results/slurm_logs/{rule}/{rule}-%j-{wildcards}.err
    #    --parsable
    #cluster-generic-status-cmd: status-sacct-robust.sh
    #cluster-generic-cancel-cmd: scancel
    #cluster-generic-cancel-nargs: 400
    #default-resources:
    #  - time="12:00:00"
    #  - mem_per_cpu=3200
    #  - tmpdir="/tmp"
    #restart-times: 2
    #max-jobs-per-second: 10
    #max-status-checks-per-second: 2
    #local-cores: 1
    #latency-wait: 60
    #cores: 800
    #jobs: 500
    #keep-going: True
    #rerun-incomplete: True

@click.command(panel = "HPC Configurations")
def hpc_lsf():
    """
    Create a template config for LSF
    
    This command creates a configuration for the LSF HPC scheduler. Writes to `stdout`.
    You will also need to install `snakemake-executor-plugin-lsf` for the HPC job submission to work.
    """
    fetch_template("hpc-lsf.yaml")
    package_absent("snakemake-executor-plugin-lsf")

@click.command(panel = "HPC Configurations")
def hpc_slurm():
    """
    Create a template config for SLURM
    
    This command creates a configuration for the SLURM HPC scheduler. Writes to `stdout`.
    You will also need to install `snakemake-executor-plugin-slurm` for the HPC job submission to work.
    """
    fetch_template("hpc-slurm.yaml")
    package_absent("snakemake-executor-plugin-slurm")

@click.command(panel = "HPC Configurations")
def hpc_googlebatch():
    """
    Create a template config for Google Batch
    
    This command creates a configuration for the Google Batch scheduler. Writes to `stdout`.
    You will also need to install `snakemake-executor-plugin-googlebatch` for the HPC job submission to work.
    """
    fetch_template("hpc-googlebatch.yaml")
    package_absent("snakemake-executor-plugin-googlebatch")


@click.command(panel =  "Input Files", context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/impute/#parameter-file")
def impute():
    """
    Create a template imputation parameter file

    With this command you can create a template parameter
    file necessary for imputation via `harpy impute`. The resulting
    file will have generic values and should be modified to be appropriate
    for your study system. Writes to `stdout`.
    """
    fetch_template("impute.tsv")
    print_notice("Modify the model parameters as needed, but [yellow bold]do not add/remove columns.")

template.add_command(impute)
template.add_command(groupings)
template.add_command(hpc_generic)
template.add_command(hpc_googlebatch)
template.add_command(hpc_lsf)
template.add_command(hpc_slurm)
