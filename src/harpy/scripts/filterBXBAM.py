import pysam
import re
import argparse

parser = argparse.ArgumentParser(prog = 'filterBXBAM.py',
                    description = 'Remove alignments from a BAM file that have a least one invalid beadtag barcode.')
parser.add_argument('-i', '--input', help = "Input bam/sam file. A corresponding index file should be in the same directory.")
parser.add_argument('-r', '--remove', type = str, help = "Barcode name or file listing BX barcodes to remove. If file, must be one BX per line.")
parser.add_argument('-v', '--valid', action='store_true', help = 'Remove alignments from a BAM file that have a least one invalid beadtag barcode.')
args = parser.parse_args()
alnfile = pysam.AlignmentFile(args.input)

if args.valid and args.remove is None:
    outfile = pysam.AlignmentFile(args.input[0:-4] + ".bx.valid.bam", "wb", template = alnfile)
    for read in alnfile.fetch():
        try:
            bx = read.get_tag("BX")
        except:
            # There is no bx tag
            continue
        # do a regex search to find X00 pattern in the BX
        if re.search("[A-Z]0{2,4}", bx):
            # if found, invalid and skipped
            continue
        else:
            outfile.write(read)

    alnfile.close()
    outfile.close()
    exit(0)

elif args.valid == False and args.remove is not None:
    try:
        with open(args.remove) as file:
            rmBX = [line.rstrip() for line in file] 
            print("USING FILE")
    except:
        rmBX = [args.remove]
    outfile = pysam.AlignmentFile(args.input[0:-4] + ".bx.filtered.bam", "wb", template = alnfile)
    for read in alnfile.fetch():
        try:
            bx = read.get_tag("BX")
        except:
            # There is no bx tag
            outfile.write(read)
            continue
        if bx in rmBX:
            # if found, skipped
            continue
        else:
            outfile.write(read)
    alnfile.close()
    outfile.close()
    exit(0)

else:
    try:
        with open(args.remove) as file:
            rmBX = [line.rstrip() for line in file] 
    except:
        rmBX = [args.remove]
    outfile = pysam.AlignmentFile(args.input[0:-4] + ".bx.valid.filtered.bam", "wb", template = alnfile)
    for read in alnfile.fetch():
        try:
            bx = read.get_tag("BX")
        except:
            # There is no bx tag
            continue
        # do a regex search to find X00 pattern in the BX
        if bx in rmBX or re.search("[A-Z]0{2,4}", bx):
            # if found, skipped
            continue
        else:
            outfile.write(read)
    alnfile.close()
    outfile.close()
    exit(0)
