import os
from shutil import rmtree
from time import sleep
import re
import rich_click as click
from rich.live import Live
from rich.panel import Panel
import subprocess
from harpy.common.printing import print_error, CONSOLE
from harpy.report.render import ReportRender

@click.command(context_settings={"help_option_names" : ['--help']}, epilog = "Documentation: https://pdimens.github.io/harpy/reports/")
@click.option('-d', '--debug', is_flag = True, help = 'Dump all of jupyterbook\'s output to the terminal')
@click.option('-h', '--headless', is_flag = True, help = 'Run the server in headless mode, with only the content server started')
@click.option('-c', '--clear-cache', is_flag = True, default = False, help = 'Remove `_build` directory prior to server launch')
@click.option('-p', '--port', type = int, help = 'Run the application server from the specified port number')
@click.option('-r', '--refresh', type = click.IntRange(min = 0, max_open=True), show_default = True, default = 0, help = 'Refresh interval, in seconds (disabled with `0`)')
@click.option('-s', '--server-port', type = int, help = 'Run the content server from the specified port number')
@click.argument('directory', required=False, type = click.Path(exists = True, file_okay = False, readable = True), nargs = 1)
def report(directory, debug, headless, clear_cache, port, server_port, refresh):
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

    # clear out the existing build dir, if present
    if os.path.isdir("_build") and clear_cache:
        rmtree("_build", ignore_errors=True)

    tracker = ReportRender(directory if directory else "")
    tracker.scan_for_reports()
    tracker.update_yaml()
    URL = ""
    myst_error = ""
    if debug:
        os.system(" ".join(cmd))
        return
    try:
        start_text = "Starting the MyST live-server[dim]…[/]" if not clear_cache else "Fetching site template[dim]…[/]"
        panel = Panel(start_text, border_style = "medium_purple4", title = "[default bold]Harpy report[/]", subtitle= "[default]Terminate server with[/] [bold yellow]ctrl+c[/]")
        with subprocess.Popen(cmd, cwd = directory, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True) as serve, Live(panel, console = CONSOLE, auto_refresh = False, transient = True) as live:
            while not URL:
                if serve.poll():
                    myst_error += "\n".join(serve.stderr.readlines())
                    raise ValueError
                _myst_output = serve.stdout.readline()
                if "Installing web libraries" in _myst_output:
                    panel.renderable = "Installing web libraries for site[dim]…[/]"
                    live.refresh()
                if "Installed web libraries" in _myst_output:
                    panel.renderable = "Starting the MyST live-server[dim]…[/]"
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
                    tracker.scan_for_reports()
                    tracker.update_yaml()
                    sleep(refresh)
                else:
                    sleep(9999)
                live.refresh()
    except KeyboardInterrupt:
        # clear the top part of the panel
        for _ in range(1):
            CONSOLE.file.write("\033[F\033[K")
        CONSOLE.file.flush()
    except ValueError:
        print_error(
            "MyST server error",
            f"The [blue]myst start[/] command exited and reported this error:\n[yellow]{myst_error.strip()}[/]"
        )
