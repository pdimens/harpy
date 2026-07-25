import os
import re
import subprocess
import sys
from shutil import rmtree
from time import sleep

import rich_click as click
from rich.live import Live
from rich.panel import Panel

from harpy.common.printing import HarpyPrint
from harpy.report.render import ReportRender
from harpy.report.static import ReportStatic
from harpy.report.utilities import check_tool
from harpy.common.cli_filetypes import IPYNBfile

@click.group(options_metavar='')
@click.command_panel("Commands", panel_styles={"border_style": "blue"})
@click.help_option('--help', hidden = True)
def report():
    """
    Render ipynb reports

    Provide an additional subcommand `live` or `static` to get more information on those
    two options. 
    """

@click.command(epilog = "Documentation: https://pdimens.github.io/harpy/reports/")
@click.option('-d', '--debug', is_flag = True, help = 'Dump all of jupyterbook\'s output to the terminal')
@click.option('-h', '--headless', is_flag = True, help = 'Run the server in headless mode, with only the content server started')
@click.option('-m', '--md', is_flag = True, help = 'Also scan for markdown files (.md)')
@click.option('-c', '--clear-cache', is_flag = True, default = False, help = 'Remove `_build` directory prior to server launch')
@click.option('-p', '--port', type = int, help = 'Run the application server from the specified port number')
@click.option('-r', '--refresh', type = click.IntRange(min = 0, max_open=True), show_default = True, default = 0, help = 'Refresh interval, in seconds (disabled with `0`)')
@click.option('-s', '--server-port', type = int, help = 'Run the content server from the specified port number')
@click.help_option('--help', hidden = True)
@click.argument('directory', required=False, type = click.Path(exists = True, file_okay = False, readable = True))
def live(directory, debug, md, headless, clear_cache, port, server_port, refresh):
    """
    Render ipynb reports as a local website

    Using MyST, all the `.ipynb` reports within Harpy-generated
    directories will be aggregated and rendered into a locally-served
    website for you to review them from a single access point. This command
    is expected to be executed within a git version-controlled directory, where
    Harpy can identify the root directory of the project, otherwise provide the
    path to a directory for Harpy to recursively scan the `.ipynb` reports.
    """
    cmd = ["jupyter", "book", "start"]
    if headless:
        cmd.append("--headless")
    if port:
        cmd += ["--port", f"{port}"]
    if server_port:
        cmd += ["--server-port", f"{server_port}"]

    hp = HarpyPrint()

    # clear out the existing build dir, if present
    if os.path.isdir("_build") and clear_cache:
        rmtree("_build", ignore_errors=True)

    tracker = ReportRender(directory if directory else os.getcwd(), md)
    tracker.scan()
    tracker.update_yaml()
    URL = ""
    myst_error = ""
    if debug:
        subprocess.run(cmd, cwd = directory, stdout = sys.stdout, stderr = sys.stderr)
        return
    try:
        start_text = "Starting the MyST live-server[dim]…[/]" if not clear_cache else "Fetching site template[dim]…[/]"
        panel = Panel(start_text, border_style = "medium_purple4", title = "[default bold]Harpy report[/]", subtitle= "[default]Terminate server with[/] [bold yellow]ctrl+c[/]")
        with (
            subprocess.Popen(cmd, cwd = directory, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True) as serve,
            Live(panel, console = hp.console, auto_refresh = False, transient = True) as live
        ):
            while not URL:
                if serve.poll():
                    myst_error += "\n".join(serve.stderr.readlines())
                    raise ValueError
                _myst_output = serve.stdout.readline()
                if "Installing web libraries" in _myst_output:
                    panel.renderable = "Installing web libraries for site[dim]…[/]"
                    live.refresh()
                if "Installed web libraries" in _myst_output:
                    panel.renderable = "Starting the report live-server[dim]…[/]"
                    live.refresh()
                _url = re.findall(r'http://\S+\s', _myst_output)
                if _url:
                    URL += _url[0].strip()
            panel.renderable = f"Report live-server running: [blue bold]{URL}[/]"
            live.refresh()
            while True:
                if serve.poll():
                    raise ValueError
                if refresh > 0:
                    tracker.scan()
                    tracker.update_yaml()
                    sleep(refresh)
                else:
                    sleep(9999)
                live.refresh()
    except KeyboardInterrupt:
        # clear the top part of the panel
        for _ in range(1):
            hp.file.write("\033[F\033[K")
        hp.file.flush()
    except ValueError:
        hp.error(
            "MyST server error",
            f"The [blue]jupyter book start[/] command exited and reported this error:\n[yellow]{myst_error.strip()}[/]"
        )


@click.command(no_args_is_help = True, context_settings={"allow_interspersed_args" : False}, epilog = "Documentation: https://pdimens.github.io/harpy/reports/")
@click.option('-d', '--debug', is_flag = True, default = False, help = 'Log process information to stderr')
@click.option('-s', '--self-contained', is_flag = True, default = False, help = 'Store all JS and CSS within the output file')
@click.help_option('--help', hidden = True)
@click.argument('notebooks', required=True, type=IPYNBfile(), nargs=-1)
def static(notebooks, debug, self_contained):
    """
    Convert ipynb reports to standalone HTML files

    Converts pre-executed .ipynb files generated by Harpy workflows into their own HTML files that do
    not require binding into a website or a live-server to be viewed. While not necessary in most cases,
    using `--self-contained` will bundle as much Javascript and CSS libraries as possible into the HTML
    file so it does not need internet access to retrieve those libraries. The resulting HTML files are
    created in the same directories as their source `.ipynb` files, but will lack the nicer features and
    formatting of a proper MyST-MD website.
    """
    check_tool(
        "jupyter",
        "It can be installed with using one of these methods:\npip: [green]pip install -U nbconvert[/]\n"
        "conda: [green]conda install -c conda-forge nbconvert[/]\n"
        "pixi: [green]pixi add nbconvert[/]"
    )
    if self_contained:
        check_tool(
            "monolith",
            "Monolith is required to flatten an HTML notebook and bundle the Javascript and CSS elements within it. "
            "Harpy does not provide it, but it can be installed using "
            "cargo ([green]cargo install monolith[/]) or by downloading and adding "
            "a pre-built binary from [blue]https://github.com/Y2Z/monolith/releases[/] to your PATH."
        )
    all_notebooks = [nb for group in notebooks for nb in group]
    n = len(all_notebooks)
    if n > 1 :
        print(f"Converting {n} notebooks into HTML files.", file = sys.stderr)
    rs = ReportStatic("", quiet = not debug, static = self_contained)
    for nb in all_notebooks:
        rs.notebook = nb
        rs.convert()

report.add_command(live)
report.add_command(static)