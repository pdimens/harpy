"""Generate shell completion scripts"""

import importlib

import rich_click as click
from click.shell_completion import get_completion_class

SHELLS = ["bash", "zsh", "fish"]
# the commands that have entry points, and where to find the ones other than `harpy` itself
PROGRAMS = {
    "harpy": None,
    "harpy-utils": "harpy.utils:cli",
    "hv": "harpy.commands.view:view",
}

@click.command(hidden = True, no_args_is_help = True)
@click.help_option('--help', hidden = True)
@click.argument('shell', required = True, type = click.Choice(SHELLS))
@click.argument('program', required = False, default = "harpy", type = click.Choice(list(PROGRAMS)))
def completion(shell, program):
    """
    Print the shell completion script

    **INTERNAL USE ONLY**. Prints the tab-completion script for `harpy` (or `harpy-utils`, or `hv`) for the
    given shell. This is run during the conda/pixi build to write the scripts to `share/`, where they are
    picked up when the environment is activated. Users should not need to run this.
    """
    if PROGRAMS[program] is None:
        # find the root group instead of importing harpy.__main__, which works
        # regardless of how harpy was invoked (entry point vs `python -m harpy`)
        command = click.get_current_context().find_root().command
    else:
        module_path, attribute = PROGRAMS[program].split(":")
        command = getattr(importlib.import_module(module_path), attribute)
    # same variable click derives from the program name at completion time, e.g. _HARPY_UTILS_COMPLETE
    complete_var = f"_{program}_COMPLETE".replace("-", "_").upper()
    click.echo(get_completion_class(shell)(command, {}, program, complete_var).source())
