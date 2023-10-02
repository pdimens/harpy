#! /usr/bin/env python3

import re
import os
import sys
import pysam
import argparse

parser = argparse.ArgumentParser(
    prog = 'assignMI.py',
    description = 
    """
    Assign an MI:i: (Molecular Identifier) tag to each barcoded 
    record based on a molecular distance cutoff. Unmapped records
    are discarded in the output. Records without a BX:Z: tag or
    with an invalid barcode (00 as one of its segments) are presevered
    but are not assigned an MI:i tag. Input file MUST BE COORDINATE SORTED.
    """,
    usage = "assignMI.py -c cutoff -i input.bam -o output.bam",
    exit_on_error = False
    )
parser.add_argument('-c','--cutoff', type=int, default = 100000, help = "Distance in base pairs at which alignments with the same barcode should be considered different molecules. (default: 100000)")
parser.add_argument('-i', '--input', help = "Input coordinate-sorted bam/sam file. If bam, a matching index file should be in the same directory.")
parser.add_argument('-o', '--output', help = "Output bam file. Will also create an index file.")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

def write_validbx(bam, alnrecord, molID):
    '''
    bam: the output bam
    alnrecord: the pysam alignment record
    molID: the "mol_id" entry of a barcode dictionary
    Formats an alignment record to include the MI tag
    and the BX at the end and writes it to the output
    bam file. Replaces existing MI tag, if exists.
    '''
    # get all the tags except MI b/c it's being replaced (if exists)
    # also remove DI because it's not necessary
    tags = [j for i,j in enumerate(alnrecord.get_tags()) if j[0] not in ['MI', 'DI']]
    # add the MI tag
    tags.append(("MI", molID))
    # find which tag index is the BX tag
    BX_idx = [i for i,j in enumerate(tags) if j[0] == 'BX'][0]
    # get the list of indices for the updated tag list
    idx = [i for i in range(len(tags))]
    # swap the BX and MI tags to make sure BX is at the end
    # b/c LEVIATHAN insists it's at the end
    tags[BX_idx], tags[-1] = tags[-1], tags[BX_idx]
    # update the record's tags
    alnrecord.set_tags(tags)
    # write record to output file
    bam.write(alnrecord)

def write_invalidbx(bam, alnrecord):
    '''
    bam: the output bam
    alnrecord: the pysam alignment record
    Formats an alignment record to include the BX 
    at the end and writes it to the output
    bam file. Removes existing MI tag, if exists.
    '''
    # get all the tags except MI b/c it's being replaced (if exists)
    # this won't write a new MI, but keeping an existing one
    # may create incorrect molecule associations by chance
    # also remove DI because it's not necessary
    tags = [j for i,j in enumerate(alnrecord.get_tags()) if j[0] not in ['MI', 'DI']]
    # find which tag index is the BX tag
    BX_idx = [i for i,j in enumerate(tags) if j[0] == 'BX'][0]
    # get the list of indices for the tag list
    idx = [i for i in range(len(tags))]
    # if BX isn't already last, make sure BX is at the end
    if tags[-1][0] != 'BX':
        # swap it with whatever is last
        tags[BX_idx], tags[-1] = tags[-1], tags[BX_idx]
        # update the record's tags
        alnrecord.set_tags(tags)
    # write record to output file
    bam.write(alnrecord)

def write_missingbx(bam, alnrecord):
    '''
    bam: the output bam
    alnrecord: the pysam alignment record
    Formats an alignment record to include invalid BX 
    at the end and writes it to the output
    bam file. Removes existing MI tag, if exists.
    '''
    # get all the tags except MI b/c it's being replaced (if exists)
    # this won't write a new MI, but keeping an existing one
    # may create incorrect molecule associations by chance
    # also remove DI because it's not necessary
    tags = [j for i,j in enumerate(alnrecord.get_tags()) if j[0] not in ['MI', 'DI']]
    tags.append(("BX", "A00C00B00D00"))
    alnrecord.set_tags(tags)
    # write record to output file
    bam.write(alnrecord)

args = parser.parse_args()

# initialize the dict
d = dict()
# chromlast keeps track of the last chromosome so we can
# clear the dict when it's a new contig/chromosome
chromlast = False
# MI is the name of the current molecule. Arbitrarily starting at 10000
MI = 9999

if os.path.exists(args.input) and args.input.lower().endswith(".sam"):
    alnfile = pysam.AlignmentFile(args.input)
elif os.path.exists(args.input) and args.input.lower().endswith(".bam"):
    if os.path.exists(args.input + ".bai"):
        alnfile = pysam.AlignmentFile(args.input)
    else:
        print(f"Error: {args.input} requires a matching {args.input}.bai index file, but one wasn\'t found.", file = sys.stderr)
        exit(1)
else:
    print(f"Error: {args.input} not found", file = sys.stderr)
    exit(1)

# iniitalize output file
#alnfile = pysam.AlignmentFile("/home/pdimens/Documents/harpy/test/bam/sample1.bam")
outfile = pysam.AlignmentFile(args.output, "wb", template = alnfile)
#outfile = pysam.AlignmentFile("/home/pdimens/Documents/harpy/test/bam/test.bam", "w", template = alnfile)

for record in alnfile.fetch():
    chrm = record.reference_name
    bp   = record.query_alignment_length
    # check if the current chromosome is different from the previous one
    # if so, empty the dict (a consideration for RAM usage)
    if chromlast != False and chrm != chromlast:
        d = dict()
    if record.is_unmapped:
        # skip, don't output
        chromlast = chrm
        continue

    try:
        bx = record.get_tag("BX")
        # do a regex search to find X00 pattern in the BX
        if re.search("[ABCD]0{2,4}", bx):
            # if found, invalid
            write_invalidbx(outfile, record)
            chromlast = chrm
            continue
    except:
        # There is no bx tag
        write_missingbx(outfile, record)
        chromlast = chrm
        continue
    
    aln = record.get_blocks()
    if not aln:
        # unaligned, skip and don't output
        chromlast = chrm
        continue

    # logic to accommodate split records 
    # start position of first alignment
    pos_start = aln[0][0]
    # end position of last alignment
    pos_end   = aln[-1][1]

    # create bx entry if it's not present
    if bx not in d.keys():
        # increment MI b/c it's a new molecule
        MI += 1 
        d[bx] = {
            "lastpos" : pos_end,
            "current_suffix": 0,
            "mol_id": MI
        }
        # write and move on
        write_validbx(outfile, record, d[bx]["mol_id"])
        chromlast = chrm
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
    if dist > args.cutoff:
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
        write_validbx(outfile, record, d[bx]["mol_id"])
        chromlast = chrm
        continue 

    if record.is_reverse or (record.is_forward and not record.is_paired):
        # set the last position to be the end of current alignment
        d[bx]["lastpos"] = pos_end

    # if it hasn't moved on by now, it's a record for an
    # existing barcode/molecule. Write the record.
    write_validbx(outfile, record, d[bx]["mol_id"])

    # update the chromosome tracker
    chromlast = chrm

alnfile.close()
outfile.close()

# index the output file
pysam.index(args.output)

# exit gracefully
exit(0)