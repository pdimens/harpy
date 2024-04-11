import argparse
import sys

parser = argparse.ArgumentParser(
    prog = 'inferSV.py',
    description = 'Create column in NAIBR bedpe output inferring the SV type from the orientation. Removes variants with FAIL flags, use optional -f argument to output FAIL variants to a separate file.',
    usage = "inferSV.py file.bedpe [-f fail.bedpe] > outfile.bedpe",
    exit_on_error = False
    )

parser.add_argument("bedfile", help = "Input bedpe file containing the output of NAIBR.")
parser.add_argument("-f", "--fail", dest = "failfile", required = False, type=str, metavar = "<failed.bedpe>", help="output variants who fail filtering into separate file")

# parser.add_argument("outfile", help = "Output bam file. A corresponding index file will be created for it.")

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

if args.failfile is not None:
    failout = open(args.failfile, "w")
    with open(args.bedfile, "r") as f:
        # first line, the header
        line = f.readline().strip().split("\t")
        line_corrected = [i.title().replace(" ", "") for i in line]
        line_corrected.append("SV")
        headerline = "\t".join(line_corrected) + "\n"
        sys.stdout.write(headerline)
        failout.write(headerline)
        # the reamining lines
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip().split("\t")
            inference = conversions[line[6]]
            line.append(inference)
            newrow = "\t".join(line) + "\n"
            if "FAIL" not in line:
                sys.stdout.write(newrow)
            else:
                failout.write(newrow)
    failout.close()
else:
    with open(args.bedfile, "r") as f:
        # first line, the header
        line = f.readline().strip().split("\t")
        line_corrected = [i.title().replace(" ", "") for i in line]
        line_corrected.append("SV")
        headerline = "\t".join(line_corrected) + "\n"
        sys.stdout.write(headerline)
        # the reamining lines
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip().split("\t")
            if "FAIL" not in line:
                inference = conversions[line[6]]
                line.append(inference)
                newrow = "\t".join(line) + "\n"
                sys.stdout.write(newrow)
