"""Harpy demultiplex workflows"""

from itertools import product
import os
from random import sample
import re
import sys
import pysam
import subprocess
import rich_click as click
from harpy.common.cli_filetypes import HPCProfile, FASTQfile
from harpy.common.cli_types_generic import SnakemakeParams
from harpy.common.misc import container_ok, safe_read
from harpy.common.printing import print_error, print_solution_with_culprits, workflow_info
from harpy.common.validations import validate_demuxschema
from harpy.common.workflow import Workflow

@click.group(options_metavar='', context_settings={"help_option_names" : ["-h", "--help"]})
def demultiplex():
    """
    Demultiplex haplotagged FASTQ files

    Check that you are using the correct haplotagging method/technology, since the different
    barcoding approaches have very different demultiplexing strategies.

    **Haplotagging Technologies**
    - `meier2021`: the original haplotagging barcode strategy
      - Meier _et al._ (2021) doi: 10.1073/pnas.2015005118
    
    """
#    **Other**
#    - `ncbi`: restore barcodes from NCBI-downloaded linked-read sequences
#    """

docstring = {
    "harpy demultiplex meier2021": [
        {
            "name": "Parameters",
            "options": ["--keep-unknown-barcodes", "--keep-unknown-samples", "--qx-rx","--schema"],
            "panel_styles": {"border_style": "blue"}
        },
        {
            "name": "Workflow Options",
            "options": ["--container", "--hpc", "--output-dir", "--quiet", "--skip-reports", "--snakemake", "--threads", "--help"],
            "panel_styles": {"border_style": "dim"}
        },
    ]
}

@click.command(deprecated = True)
def gen1():
    """
    Use `meier2021` instead
    """
    print("This workflow has been renamed \"meier2021\"-- please use that instead. This warning will be removed in the next minor Harpy version and will only return an error.")
    sys.exit(1)

@click.command(no_args_is_help = True, context_settings=dict(allow_interspersed_args=False), epilog = "Documentation: https://pdimens.github.io/harpy/workflows/demultiplex/")
@click.option('-u', '--keep-unknown-samples',  is_flag = True, default = False, help = 'Keep a separate file of reads with recognized barcodes but don\'t match any sample in the schema')
@click.option('-b', '--keep-unknown-barcodes',  is_flag = True, default = False, help = 'Keep a separate file of reads with unrecognized barcodes')
@click.option('-q', '--qx-rx', is_flag = True, default = False, help = 'Include the `QX:Z` and `RX:Z` tags in the read header')
@click.option('-t', '--threads', default = 4, show_default = True, type = click.IntRange(2,999, clamp = True), help = 'Number of threads to use')
@click.option('-o', '--output-dir', type = click.Path(exists = False, resolve_path = True), default = "Demultiplex", show_default=True,  help = 'Output directory name')
@click.option('--container',  is_flag = True, default = False, help = 'Use a container instead of conda', callback=container_ok)
@click.option('--setup-only',  is_flag = True, hidden = True, default = False,  help = 'Setup the workflow and exit')
@click.option('--hpc',  type = HPCProfile(), help = 'HPC submission YAML configuration file')
@click.option('--quiet', show_default = True, default = 0, type = click.Choice([0, 1, 2]), help = '`0` all output, `1` show one progress bar, `2` no output')
@click.option('--skip-reports',  is_flag = True, show_default = True, default = False, help = 'Don\'t generate HTML reports')
@click.option('--snakemake', type = SnakemakeParams(), help = 'Additional Snakemake parameters, in quotes')
@click.argument('schema', required = True, type=click.Path(exists=True, dir_okay=False, readable=True, resolve_path=True))
@click.argument('R12_FQ', required=True, type=FASTQfile(dir_ok= False), nargs=2)
@click.argument('I12_FQ', required=True, type=FASTQfile(dir_ok= False), nargs=2)
def meier2021(r12_fq, i12_fq, output_dir, schema, qx_rx, keep_unknown_samples, keep_unknown_barcodes, threads, snakemake, skip_reports, quiet, hpc, container, setup_only):
    """
    Demultiplex FASTQ files haplotagged with the Meier _et al._ 2021 protocol

    Use the R1, R2, I2, and I2 FASTQ files provided by the sequencing facility as inputs after the options and schema (4 files, in that exact order). 
    The `SCHEMA` must have **no header** (i.e. no column names) and be in the format of `sample`\\<TAB\\>`barcode`,
    where `barcode` is the barcode segment associated with the sample ID (.e.g. `C01`, `C02`, etc.). Use `--qx-rx` to add the 
    `QX:Z` (barcode PHRED scores) and `RX:Z` (nucleotide barcode) tags in the sequence headers. These tags aren't used by any
    subsequent analyses, but may be useful for your own diagnostics. 
    """
    workflow = Workflow("demultiplex_meier2021", "demultiplex_meier2021.smk", output_dir, quiet) 
    workflow.setup_snakemake(container, threads, hpc, snakemake)
    workflow.conda = ["demultiplex", "qc"]
    
    ## checks and validations ##
    multi_id = validate_demuxschema(schema)

    workflow.config = {
        "workflow" : workflow.name,
        "retain" : {
            "qx_rx" : qx_rx,
            "barcodes" : keep_unknown_barcodes,
            "samples" : keep_unknown_samples,
        },
        "snakemake" : {
            "log" : workflow.snakemake_log,
            "absolute": workflow.snakemake_cmd_absolute,
            "relative": workflow.snakemake_cmd_relative,
        },
        "conda_environments" : workflow.conda,
        "reports" : {
            "skip": skip_reports
        },
        "inputs" : {
            "demultiplex_schema" : schema,
            "R1": r12_fq[0],
            "R2": r12_fq[1],
            "I1": i12_fq[0],
            "I2": i12_fq[1]
        }
    }
    
    workflow.start_text = workflow_info(
        ("Barcode Design:", "Meier [italic]et al.[/] 2021"),
        ("Demultiplex Schema:", os.path.basename(schema)),
        ("Multi-ID Samples:", "[bold yellow]Yes[/] [dim](assuming this was intentional)[/]" if multi_id else "No"),
        ("Include QX/RX tags", "Yes" if qx_rx else "No"),
        ("Output Folder:", os.path.basename(output_dir) + "/")
    )

    workflow.initialize(setup_only)

@click.command(hidden = True, no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/harpy/ncbi")
@click.option('-m', '--barcode-map', type=click.Path(exists=True, readable=True, resolve_path=True), help = 'Map of nucleotide-to-barcode conversion')
@click.option('-s', '--barcode-style', type=click.Choice(["haplotagging","nucleotide","tellseq","10x", "stlfr"], case_sensitive = False), help = 'Style format for barcodes in output')
@click.option('-p', '--prefix', required=True, type = str, help = "Output file name prefix")
@click.argument('r1_fq', required=True, type=click.Path(exists=True, readable=True, resolve_path=True), nargs=1)
@click.argument('r2_fq', required=True, type=click.Path(exists=True, readable=True, resolve_path=True), nargs=1)
def ncbi(prefix, r1_fq, r2_fq, barcode_map, barcode_style):
    """
    Restore linked-read barcodes to the sequence header

    This demultiplexing strategy is the complement to `harpy convert ncbi`, restoring barcodes to the
    sequence headers and removing barcodes and spacers from the sequences. By default, it will keep
    barcodes in nucleotide format unless a `--barcode-map`/`-m` is provided with specific conversions in `nucleotide<tab>barcode` format
    or a `--barcode-style`/`-s` is given. The output files will have barcodes encoded in the Standard
    configuration, i.e. as a `BX:Z` tag along with a `VX:i` tag indicating barcode validation. Requires
    a `PREFIX` to name the output files.

    | --barcode-style              | format                | example             |
    |:-----------------------------|:----------------------|:--------------------|
    | `nucleotide`/`tellseq`/`10x` | nucleotides (default) | `BX:Z:ATAGGACGAAGA` |
    | `haplotagging`               | AxxCxxBxxDxx          | `BX:Z:A01C93B56D11` |
    | `stlfr`                      | 1_2_3                 | `BX:Z:154_211_934`  |
    """
    if barcode_map and barcode_style:
        print_error("invalid options", "[blue]--barcode-map[/] and [blue]--barcode-style[/] cannot be used together.")
        sys.exit(1)
    conv_dict = {}
    if barcode_map:
        with safe_read(barcode_map) as bcmap:
            for i,j in enumerate(bcmap,1):
                try:
                    nuc,bx = j.split()
                    if not re.search(r"^[ATCGN]+$", nuc):
                        print_error("Bad file format", f"The file provided to [blue]--barcode-map[/] requires nucleotide barcodes in the first column, but characters other than [green]ATCGN[/] were found in row [bold]{i}[/]")
                        print_solution_with_culprits("Make sure the mapping file you are providing is in the format:\n[green]nucleotides[/][dim]<tab or space>[/][green]new_barcode[/]", f"Contents of row {i}")
                        click.echo(j.strip())
                        sys.exit(1)
                except ValueError:
                    print_error("Bad file format", f"The file provided to [blue]--barcode-map[/] expects two entries per row separated by a whitespace, but a different amount was found in row [bold]{i}[/]")
                    print_solution_with_culprits("Make sure the mapping file you are providing is in the format:\n[green]nucleotides[/][dim]<tab or space>[/][green]new_barcode[/]", f"Contents of row {i}")
                    click.echo(j.strip())
                    sys.exit(1)
                conv_dict[nuc] = bx

    if barcode_style == "stlfr":
        bc_generator = product(*[sample(range(1,1538), 1537) for i in range(3)])
        INVALID = "0_0_0"
        SEP = "_"
    elif barcode_style == "haplotagging":
        bc_generator = product(
            ["A" + str(i).zfill(2) for i in sample(range(1,97), 96)],
            ["C" + str(i).zfill(2) for i in sample(range(1,97), 96)],
            ["B" + str(i).zfill(2) for i in sample(range(1,97), 96)],
            ["D" + str(i).zfill(2) for i in sample(range(1,97), 96)]
        )
        INVALID = "A00C00B00D00"
        SEP = ""
    else:
        INVALID = "N"*18
        SEP = ""

    def format_bc(bc):
        return SEP.join(str(i) for i in bc)

    def bx_position(rec):
        # search first 30 bases and return the INDEX of where the NNNNN spacer ends
        _bx = re.search(r"[ATCGN]*NNNNN", rec.sequence[:30])
        if not _bx:
            return None
        else:
            return _bx.end()

    def is_barcode(string):
        '''
        Returns True if the quality string is all I and ends with !!!!!, which is the expected encoding
        Also considers just !!!!! without any prefix as valid, i.e. it's just a spacer with no barcode
        '''
        if len(string) < 5:
            return False
        if not string.endswith("!"*5):
            return False
        if string == "!"*5:
            return True
        return all(c == "I" for c in string[:-5])

    sys.stderr.write("\t".join(["File", "Reads_Demultiplexed", "Reads_Ambiguous"]) + "\n")
    for i,fq in enumerate([r1_fq, r2_fq],1):
        sys.stderr.write(os.path.basename(fq))
        with (
            pysam.FastqFile(fq, "r") as in_fq,
            open(f"{prefix}.R{i}.fq.gz", "wb") as out_fq,
            open(f"{prefix}.ambiguous.R{i}.fq.gz", "wb") as ambig_fq
        ):
            gzip = subprocess.Popen(["gzip"], stdin = subprocess.PIPE, stdout = out_fq)
            gzip_ambig = subprocess.Popen(["gzip"], stdin = subprocess.PIPE, stdout = ambig_fq)
            
            AMBIG_TOTAL = 0
            DEMUX_TOTAL = 0

            for record in in_fq:
                AMBIGUOUS = False
                _bx_idx = bx_position(record)
                if _bx_idx:
                    _bx = record.sequence[:_bx_idx]
                    _qual = record.quality[:_bx_idx]
                    if not is_barcode(_qual):
                        # the format is different that all I followed by 5 !
                        AMBIGUOUS = True
                    else:
                        record.sequence = record.sequence[_bx_idx:]
                        record.quality = record.quality[_bx_idx:]
                        _bx = _bx.removesuffix("N"*5)
                        # or not _bx safeguards against just a spacer
                        invalid = "N" in _bx or not _bx
                        if _bx in conv_dict:
                            new_bx = conv_dict[_bx]
                        elif not barcode_map:
                            # only create a new barcode if a map wasn't provided
                            if barcode_style not in ["stlfr", "haplotagging"]:
                                new_bx = _bx if _bx else INVALID
                            else:
                                try:
                                    new_bx = format_bc(next(bc_generator)) if not invalid else INVALID
                                except StopIteration:
                                    print_error("no more barcodes", f"There are more unique barcodes in the input files than {barcode_style} can support. Consider using [blue bold]tellseq[/] format to retain barcodes as nucleotides.")
                            if _bx:
                                conv_dict[_bx] = new_bx
                        elif barcode_map and invalid:
                            new_bx = INVALID
                        else:
                            # the barcode wasn't in the provided map nor was it invalid, retain the barcode as-is
                            AMBIGUOUS = True
                            new_bx = _bx

                        record.comment = "VX:i:0" if invalid else "VX:i:1"
                        record.comment += f"\tBX:Z:{new_bx}"
                else:
                    AMBIGUOUS = True

                if AMBIGUOUS:
                    AMBIG_TOTAL += 1
                    gzip_ambig.stdin.write(str(record).encode("utf-8") + b"\n")
                else:
                    DEMUX_TOTAL += 1
                    gzip.stdin.write(str(record).encode("utf-8") + b"\n")
            sys.stderr.write(f"\t{DEMUX_TOTAL}\t{AMBIG_TOTAL}\n")
        

demultiplex.add_command(meier2021)
demultiplex.add_command(ncbi)
demultiplex.add_command(gen1)
