#! /usr/bin/env python
"""Create column in NAIBR bedpe output inferring the SV type from the orientation"""
import argparse
import sys

parser = argparse.ArgumentParser(
    prog = 'infer_sv.py',
    description = 'Create column in NAIBR bedpe output inferring the SV type from the orientation. Removes variants with FAIL flags, use optional -f argument to output FAIL variants to a separate file.',
    usage = "infer_sv.py file.bedpe [-f fail.bedpe] > outfile.bedpe",
    exit_on_error = False
    )

parser.add_argument("bedfile", help = "Input bedpe file containing the output of NAIBR.")
parser.add_argument("-f", "--fail", dest = "failfile", type=str, metavar = "fail.bedpe", help="output variants who fail filtering into separate file")

if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

conversions = {
    "+-": "deletion",
    "--": "inversion",
    "++": "inversion",
    "-+": "duplication"
    }

if args.failfile:
    with open(args.bedfile, "r", encoding="utf-8") as f, open(args.failfile, "w", encoding="utf-8") as failout:
        # first line, the header
        line = f.readline().strip().split("\t")
        line_corrected = [i.title().replace(" ", "") for i in line]
        line_corrected.append("SV")
        HEADERLINE = "\t".join(line_corrected) + "\n"
        sys.stdout.write(HEADERLINE)
        failout.write(HEADERLINE)
        # the reamining lines
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip().split("\t")
            inference = conversions[line[6]]
            line.append(inference)
            NEWROW = "\t".join(line) + "\n"
            if "FAIL" not in line:
                sys.stdout.write(NEWROW)
            else:
                failout.write(NEWROW)
else:
    with open(args.bedfile, "r", encoding="utf-8") as f:
        # first line, the header
        line = f.readline().strip().split("\t")
        line_corrected = [i.title().replace(" ", "") for i in line]
        line_corrected.append("SV")
        HEADERLINE = "\t".join(line_corrected) + "\n"
        sys.stdout.write(HEADERLINE)
        # the reamining lines
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip().split("\t")
            if "FAIL" not in line:
                inference = conversions[line[6]]
                line.append(inference)
                NEWROW = "\t".join(line) + "\n"
                sys.stdout.write(NEWROW)
