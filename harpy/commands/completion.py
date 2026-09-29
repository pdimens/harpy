"""Generate shell completion scripts"""

import rich_click as click
from click.shell_completion import get_completion_class

SHELLS = ["bash", "zsh", "fish"]

@click.command(hidden = True, no_args_is_help = True)
@click.help_option('--help', hidden = True)
@click.argument('shell', required = True, type = click.Choice(SHELLS))
def completion(shell):
    """
    Print the shell completion script

    **INTERNAL USE ONLY**. Prints the tab-completion script for `harpy` for the given shell.
    This is run during the conda/pixi build to write the scripts to `share/harpy/`, which
    are sourced when the environment is activated. Users should not need to run this.
    """
    # find the root group instead of importing harpy.__main__, which works
    # regardless of how harpy was invoked (entry point vs `python -m harpy`)
    root = click.get_current_context().find_root().command
    completer = get_completion_class(shell)(root, {}, "harpy", "_HARPY_COMPLETE")
    click.echo(completer.source())
