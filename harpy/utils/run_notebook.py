import os
import random
import string
import sys
import tempfile
from datetime import datetime

from papermill.cli import _resolve_type

import click
import papermill as pm
from traitlets.config import Config

UID = ''.join(random.choices(string.ascii_letters + string.digits, k=15))

def _process(lines, text):
    _date = datetime.now().strftime('%Y-%m-%d')
    text = list(text)
    for line in lines:
        if line.startswith("Ctrl click to launch") or "Starting kernel" in line:
            continue
        if 'PLACEHOLDER' in line:
            if not text:
                sys.stderr.write("ERROR: more PLACEHOLDER text than replacement text provided\n")
                sys.exit(1)
            line = line.replace("PLACEHOLDER", text.pop(0))
        elif "9999-12-31" in line:
            line = line.replace("9999-12-31", _date)
        elif "injected-parameters" in line:
            line = line.replace('"injected-parameters"', '"injected-parameters",\n"remove-cell"')
        if "placeholder" in line:
            line = line.replace("placeholder", UID)
        sys.stdout.write(line)


@click.command(no_args_is_help=True)
@click.option("-k", "--kernel", required=True, type=str)
@click.option("-p", "--parameter", "params", multiple=True, nargs=2, type=str)
@click.argument("notebook", required=True, type=click.Path(exists=True, dir_okay=False))
@click.argument("text", nargs=-1, type=str)
@click.help_option('--help', hidden=True)
def run_notebook(kernel, params, notebook, text):
    """
    INTERNAL use only. Execute a notebook with papermill (IPC kernel transport) and process it.
    Writes to stdout. Everything prior to the input notebook accepts `papermill` style inputs of
    `-p <parameterName> <parameterValue>`, as many as you want. Everything after the notebook
    are the text replacements for notebook PLACEHOLDERs, in the order which PLACEHOLDER appears.

    run-notebook -k ipython-harpy -p indir path input.ipynb arg1 arg2... > output.ipynb
    """
    # ponytail: IPC sockets avoid the TCP port race between concurrent kernels on one node
    tempfile.tempdir = f".harpyreports/{UID}"
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as tmpdir, tempfile.NamedTemporaryFile("r", suffix=".ipynb") as tmp:
        os.environ.setdefault("JUPYTER_RUNTIME_DIR", os.path.join(tmpdir, "rt"))
        os.environ.setdefault("IPYTHONDIR", os.path.join(tmpdir, "ipy"))
        c = Config()
        c.KernelManager.transport = "ipc"
        # absolute per-run socket path; default "kernel-ipc" is relative to CWD (shared, Lustre)
        c.KernelManager.ip = os.path.join(tmpdir, "kernel")
        pm.execute_notebook(
            input_path = notebook,
            output_path = tmp.name,
            parameters = {k: _resolve_type(v) for k, v in params},
            kernel_name = kernel,
            start_timeout = 120,
            progress_bar = False,
            log_output = False,
            config= c ,
        )
        _process(tmp, text)