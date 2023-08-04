from .harpymisc import getnames_err
import os
import glob
import rich_click as click
import re
import sys
import subprocess

try:
    harpypath = '{CONDA_PREFIX}'.format(**os.environ) + "/bin"
except:
    pass

@click.command(no_args_is_help = True)
@click.option('-p', '--popgroup', required = False, type=click.Path(exists=True), metavar = "Input folder Path", help = 'Create generic sample-group file using existing sample file names (fq.gz or bam)')
@click.option('-s', '--stitch-params', type=str, metavar = "Output file name", help = 'Create template STITCH parameter file')
@click.option('-h', '--hpc', type = click.Choice(["slurm", "sge"], case_sensitive = False), help = 'Create HPC scheduling profile')
def extra(popgroup, stitch_params, hpc):
    """
    Create various optional/necessary input files

    With this command you can generate a sample grouping file (for variant calling),
    a templace STITCH parameter file (for imputation), and a HPC profile for running
    Harpy on a cluster. You can use any combination of options at a time. 
    """
    if popgroup is not None:
        print('\033[1m' + "<><> Sampling Grouping File <><>" + '\033[0m', file = sys.stderr)
        try:
            samplenames = getnames_err(popgroup, '.bam')
        except:
            flist = [os.path.basename(i) for i in glob.iglob(f"{popgroup}/*") if not os.path.isdir(i)]
            r = re.compile(".*\.f(?:ast)?q(?:\.gz)?$", flags=re.IGNORECASE)
            fqlist = list(filter(r.match, flist))
            bn_r = r"[\.\_][RF](?:[12])?(?:\_00[1-9])*\.f(?:ast)?q(?:\.gz)?$"
            samplenames = set([re.sub(bn_r, "", i, flags = re.IGNORECASE) for i in fqlist])
            if len(samplenames) < 1:
                print(f"\033[1;33mERROR:\033[00m No files ending in fq.gz, fastq.gz, or .bam found in {popgroup}", file = sys.stderr)
                sys.exit(1)

        print(f"Samples detected in {popgroup}: " + str(len(samplenames)), file = sys.stderr)
        fout = "samples.groups"
        if exists("samples.groups"):
            overwrite = input("File \'samples.groups\' already exists, overwrite (no|yes)?  ").lower()
            if (overwrite == "no") or (overwrite == "n"):
                fout = input("Please suggest a different name for the output file: ")
            elif (overwrite == "yes") or (overwrite == "y"):
                fout = "samples.groups"
        with open(fout, "w") as file:
            for i in samplenames:
                file.write(i + '\tpop1\n') 
        print('Created sample population grouping file: ' + fout + '\nPlease review it, as all samples have been grouped into a single population\n', file = sys.stderr)

    if stitch_params is not None:
        print('\033[1m' + "<><> STITCH Parameter File <><>" + '\033[0m', file = sys.stderr)
        with open(stitch_params, "w") as file:
            file.write('model\tuseBX\tk\ts\tnGen\npseudoHaploid\tTRUE\t10\t5\t50\npseudoHaploid\tTRUE\t10\t1\t50\npseudoHaploid\tTRUE\t15\t10\t100')
        print(f"Created example parameter file: {stitch_params}", file = sys.stderr)
        print("Modify the model parameters as needed, but " + '\033[1m' + "DO NOT" + '\033[0m' + " add/remove columns", file = sys.stderr)

    if hpc is not None:
        print('\033[1m' + "<><> HPC Profile <><>" + '\033[0m')
        subprocess.run(["hpc_profile.py", hpc])
