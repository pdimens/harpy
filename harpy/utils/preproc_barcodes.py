import click
import json as js
import pysam
import sys

@click.command(hidden = True, no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/preprocess/")
@click.argument("json", required = True, type=click.Path(exists = True, dir_okay=False, resolve_path=True))
@click.argument("bam", required = True, type=click.Path(exists = True, dir_okay=False, resolve_path=True))
@click.help_option('--help', hidden = True)
def preproc_barcodes(json, bam):
    """
    Rename and record demultiplexed linked read barcodes coming out of Pheniqs

    INTERNAL script used in harpy preprocess, not intended for typical use.
    """
    def reconstruct_barcode(nuc_bc: str):
        segments = nuc_bc.split("-")
        A = "A" + bc.get(segments[1], "00")
        B = "B" + bc.get(segments[2], "00")
        C = "C" + bc.get(segments[3], "00")
        D = "D" + stagger.get(segments[0], "00")
        valid = not any([A=="A00", B=="B00", C=="C00", D=="D00"])
        return f"{A}{C}{B}{D}", valid
        
    with open(json, 'r') as file:
        data = js.load(file)
    stagger = {}; bc = {}
    for k,v in data['decoder']['stagger']['codec'].items():
        stagger[v['barcode'][0]] = k.removeprefix("@")
    for k,v in data['decoder']['segment']['codec'].items():
        bc[v['barcode'][0]] = k.removeprefix("@")
    del data

    with (
        pysam.AlignmentFile(bam, require_index=False, check_sq = False) as bam,
        pysam.AlignmentFile(sys.stdout.buffer, 'w', template=bam, add_sam_header=False) as out
        ):
        for record in bam.fetch(until_eof=True):
            if record.has_tag("RX"):
                bx, vx = reconstruct_barcode(record.get_tag("RX"))
                record.set_tag("BX", bx, 'Z')
                record.set_tag("VX", vx, 'i')
            out.write(record)