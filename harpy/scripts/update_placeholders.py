#! /usr/bin/env python3

import argparse
import sys

def main():
    parser = argparse.ArgumentParser(
        prog='update_placeholders',
        description='Replace all instances of PLACEHOLDER in a Jupyter notebook with the input arguments, sequentially. In other words, the first instance is replaced with the first argument, second with the second, etc.',
        usage = "update_placeholders arg1 arg2... input.ipynb > output.ipynb",
        )
    parser.add_argument('text', nargs='+', help = 'text items to replace PLACEHOLDER, separated by spaces')
    parser.add_argument('notebook', type=argparse.FileType('r'), help = "Input Jupyter notebook")

    args = parser.parse_args()
    if not args.notebook:
        print("ERROR: Requires a final positional argument to be a file")
        parser.print_help(sys.stderr)
        sys.exit(1)

    for line in args.notebook:
        if 'PLACEHOLDER' in line:
            if args.text:
                line = line.replace("PLACEHOLDER", args.text.pop(0))
            else:
                print(f"ERROR: more PLACEHOLDER text than replacement text provided", file=sys.stderr)
                sys.exit(1)
        sys.stdout.write(line)