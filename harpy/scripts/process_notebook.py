#! /usr/bin/env python3

import argparse
from datetime import datetime
import sys
import random
import string

def main():
    parser = argparse.ArgumentParser(
        prog='process-notebook',
        description='Replace all instances of PLACEHOLDER in a Jupyter notebook with the input arguments, sequentially. In other words, the first instance is replaced with the first argument, second with the second, etc. Also replaces the date-time placeholder with the actual date, adds the remove-cell tag to injected paramters, and replaces lowercase instances of placeholder with a 15 digit random alphanumeric string.',
        usage = "process-notebook arg1 arg2... input.ipynb > output.ipynb",
        )
    parser.add_argument('text', nargs='+', help = 'text items to replace PLACEHOLDER, separated by spaces')
    parser.add_argument('notebook', type=argparse.FileType('r'), help = "Input Jupyter notebook")

    args = parser.parse_args()
    if not args.notebook:
        print("ERROR: Requires a final positional argument to be a file")
        parser.print_help(sys.stderr)
        sys.exit(1)

    _date = datetime.now().strftime('%Y-%m-%d')
    random_string = ''.join(random.choices(string.ascii_letters + string.digits, k=15))

    for line in args.notebook:
        if 'PLACEHOLDER' in line:
            if args.text:
                line = line.replace("PLACEHOLDER", args.text.pop(0))
            else:
                print(f"ERROR: more PLACEHOLDER text than replacement text provided", file=sys.stderr)
                sys.exit(1)
        elif "9999-12-31" in line:
            line = line.replace("9999-12-31", _date)
        elif "injected-parameters" in line:
            line = line.replace('"injected-parameters"', '"injected-parameters",\n"remove-cell"')
        if "placeholder" in line:
            line = line.replace("placeholder", random_string)
        sys.stdout.write(line)