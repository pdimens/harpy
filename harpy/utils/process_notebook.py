import random
import string
import sys
from datetime import datetime

import click


@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/preprocess/")
#@click.option("--pdf", required = True, type=str)
@click.argument("notebook", required = True, type=click.File())
@click.argument("text", nargs = -1, type = str)
@click.help_option('--help', hidden = True)
def process_notebook(notebook, text):
    """
    Replace placeholder text in jupyter notebooks

    INTERNAL use only.Replace all instances of PLACEHOLDER in a Jupyter notebook with the input arguments,
    sequentially. In other words, the first instance is replaced with the first argument, second with the
    second, etc. Also replaces the date-time placeholder with the actual date, adds the remove-cell tag to
    injected paramters, and replaces lowercase instances of placeholder with a 15 digit random alphanumeric string.
    The `--pdf` option replaces `PDF_FILE` of `PDF_FILE-placeholder.pdf`.Writes to stdout. Example:

    process-notebook --pdf aloe input.ipynb arg1 arg2... > output.ipynb
    """
    _date = datetime.now().strftime('%Y-%m-%d')
    uid = ''.join(random.choices(string.ascii_letters + string.digits, k=15))
    text = list(text)
    for line in notebook:
        if line.startswith("Ctrl click to launch"):
            continue
        if 'PLACEHOLDER' in line:
            if text:
                line = line.replace("PLACEHOLDER", text.pop(0))
            else:
                sys.stderr.write(f"ERROR: more PLACEHOLDER text than replacement text provided\n")
                sys.exit(1)
        elif "9999-12-31" in line:
            line = line.replace("9999-12-31", _date)
        elif "injected-parameters" in line:
            line = line.replace('"injected-parameters"', '"injected-parameters",\n"remove-cell"')
        if "placeholder" in line:
            line = line.replace("placeholder", uid)
        sys.stdout.write(line)
