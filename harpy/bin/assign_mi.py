#! /usr/bin/env python
"""assign molecular identifier (MI) tags to alignments based on distance and barcode"""
import re
import os
import sys
import argparse
import pysam

parser = argparse.ArgumentParser(
    prog = 'assign_mi.py',
    description =
    """
    Assign an MI:i: (Molecular Identifier) tag to each barcoded
    record based on a molecular distance cutoff. Unmapped records
    are discarded in the output. Records without a BX:Z: tag or
    with an invalid barcode (00 as one of its segments) are presevered
    but are not assigned an MI:i tag. Input file MUST BE COORDINATE SORTED.
    """,
    usage = "assign_mi.py -c cutoff input.bam > output.bam",
    exit_on_error = False
    )

parser.add_argument('-c','--cutoff', type=int, default = 100000, help = "Distance in base pairs at which alignments with the same barcode should be considered different molecules. If 0, then alignment distance is ignored. (default: 100000)")
parser.add_argument('input', help = "Input coordinate-sorted bam/sam file. If bam, a matching index file should be in the same directory.")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()
bam_input = args.input

if not os.path.exists(bam_input):
    parser.error(f"{bam_input} was not found")

# Do a quick scan of the file until the first barcode
# to assess what kind of linked read tech it is and set
# the appropriate kind of MISSING barcode
with pysam.AlignmentFile(bam_input) as alnfile:
    for record in alnfile.fetch(until_eof = True):
        try:
            bx = record.get_tag("BX")
            if re.search(r"[ATCGN]+", bx):
                # tellseq
                MISSING = "N" * len(bx)
            elif re.search(r"\d+_\d+_\d+", bx):
                # stlfr
                MISSING = "0_0_0"
            else:
                MISSING = "A00C00B00D00"
            break
        except KeyError:
            continue


def write_validbx(bam, alnrecord, mol_id):
    '''
    bam: the output bam
    alnrecord: the pysam alignment record
    mol_id: the "mol_id" entry of a barcode dictionary
    Formats an alignment record to include the MI tag
    and the BX at the end and writes it to the output
    bam file. Replaces existing MI tag, if exists.
    '''
    # get all the tags except MI b/c it's being replaced (if exists)
    # will manually parse BX, so omit that too
    # also remove DI because it's not necessary
    tags = [j for j in alnrecord.get_tags() if j[0] not in ['MI', 'DI', 'BX']]
    tags.append(("MI", mol_id))
    _bx = alnrecord.get_tag("BX")
    if "-" in _bx:
        # it's been deconvolved, set it to a DX tag
        tags.append(("DX", _bx))
    bx_clean = _bx.split("-")[0]
    tags.append(("BX", bx_clean))
    alnrecord.set_tags(tags)
    # write record to output file
    bam.write(alnrecord)

def write_invalidbx(bam, alnrecord):
    '''
    bam: the output bam
    alnrecord: the pysam alignment record
    Formats an alignment record to include the BX 
    at the end and writes it to the output
    bam file. Keeps existing MI tag if present.
    '''
    # will not keep an existing MI tag if present
    # also remove DI because it's not necessary
    tags = [j for j in alnrecord.get_tags() if j[0] not in ['MI', 'DI', 'BX']]
    _bx = alnrecord.get_tag("BX")
    # if hyphen is present, it's been deconvolved and shouldn't have been
    # and rm the hyphen part
    bx_clean = _bx.split("-")[0]
    tags.append(("BX", bx_clean))
    # update the record's tags
    alnrecord.set_tags(tags)
    # write record to output file
    bam.write(alnrecord)

def write_missingbx(bam, alnrecord, missing):
    '''
    bam: the output bam
    alnrecord: the pysam alignment record
    Formats an alignment record to include invalid BX 
    at the end and writes it to the output
    bam file. Removes existing MI tag, if exists.
    '''
    # removes MI and DI tags, writes new BX tag
    tags = [j for j in alnrecord.get_tags() if j[0] not in ['MI', 'DI', 'BX']]
    tags.append(("BX", missing))
    alnrecord.set_tags(tags)
    # write record to output file
    bam.write(alnrecord)

# initialize the dict
d = {}
# LAST_CONTIG keeps track of the last contig so we can
# clear the dict when it's a new contig/chromosome
LAST_CONTIG = False
# MI is the name of the current molecule, starting at 1 (0+1)
MI = 0

# iniitalize input/output files
with (
    pysam.AlignmentFile(bam_input) as alnfile,
    pysam.AlignmentFile(sys.stdout.buffer, "wb", template = alnfile) as outfile
):
    for record in alnfile.fetch(until_eof = True):
        chrm = record.reference_name
        bp   = record.query_alignment_length
        # check if the current chromosome is different from the previous one
        # if so, empty the dict (a consideration for RAM usage)
        if LAST_CONTIG and chrm != LAST_CONTIG:
            d = {}
        if record.is_unmapped:
            # skip, don't output
            LAST_CONTIG = chrm
            continue

        try:
            bx = record.get_tag("BX")
            if record.get_tag("VX") == 0:
                # VX:i:0 is invalid
                write_invalidbx(outfile, record)
                LAST_CONTIG = chrm
                continue
        except KeyError:
            # There is no bx tag
            write_missingbx(outfile, record, MISSING)
            LAST_CONTIG = chrm
            continue

        aln = record.get_blocks()
        if not aln:
            # unaligned, skip and don't output
            LAST_CONTIG = chrm
            continue

        # logic to accommodate split records
        # start position of first alignment
        pos_start = aln[0][0]
        # end position of last alignment
        pos_end   = aln[-1][1]

        # create bx entry if it's not present
        if bx not in d:
            # increment MI b/c it's a new molecule
            MI += 1
            d[bx] = {
                "lastpos" : pos_end,
                "current_suffix": 0,
                "mol_id": MI
            }
            # write and move on
            write_validbx(outfile, record, MI)
            LAST_CONTIG = chrm
            continue

        # store the original barcode as `orig` b/c we might need to suffix it
        orig = bx
        # if there is a suffix, append it to the barcode name
        if d[orig]["current_suffix"] > 0:
            bx = orig + "." + str(d[orig]["current_suffix"])

        # distance from last alignment = current aln start - previous aln end
        dist = pos_start - d[bx]["lastpos"]
        # if the distance between alignments is > cutoff, it's a different molecule
        # so we'll +1 the suffix of the original barcode and relabel this one as
        # BX + suffix. Since it's a new entry, we initialize it and move on
        if args.cutoff > 0 and dist > args.cutoff:
            # increment MI b/c it's a new molecule
            MI += 1
            # increment original barcode's suffix
            d[orig]["current_suffix"] += 1
            bx = orig + "." + str(d[orig]["current_suffix"])
            # add new entry for new suffixed barcode with unique MI
            d[bx] = {
                "lastpos" : pos_end,
                "current_suffix": 0,
                "mol_id": MI
            }
            # write and move on
            write_validbx(outfile, record, MI)
            LAST_CONTIG = chrm
            continue

        if record.is_reverse or (record.is_forward and not record.is_paired):
            # set the last position to be the end of current alignment
            d[bx]["lastpos"] = pos_end

        # if it hasn't moved on by now, it's a record for an
        # existing barcode/molecule. Write the record.
        write_validbx(outfile, record, d[bx]["mol_id"])

        # update the chromosome tracker
        LAST_CONTIG = chrm
