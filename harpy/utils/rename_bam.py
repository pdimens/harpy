import os
import sys
import click

@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/preprocess/")
@click.help_option('--help', hidden = True)
@click.option("-d", "--delete", is_flag = True, default = False, help="delete the original file")
@click.argument("name", required = True, type=str)
@click.argument("input", required = True, type=click.Path(exists = True, dir_okay=False, resolve_path=True))
def rename_bam(name, input, delete):
    """
    Rename a SAM/BAM file and modify the @RG tag
    
    This is the proper way to rename a SAM/BAM file to reflect the change for both ID and SM.
    This process creates a new file \'newname.bam\' and you may use -d to delete the original file. Requires samtools.
    """
    outdir = os.path.dirname(input)
    JOB_STATUS = os.system(
        f"samtools addreplacerg -r \"ID:{name}\\tSM:{name}\" -o {outdir}/{name}.bam {input}"
        )

    if JOB_STATUS != 0:
        sys.stderr.write("samtools addreplacerg failed with an error\n")  
        sys.exit(JOB_STATUS)
    else:
        if delete:
            try:
                os.remove(input)
            except OSError:
                sys.stderr.write(f"Failed to delete {input}, but otherwise samtools was successful.\n")
        sys.exit(JOB_STATUS)
