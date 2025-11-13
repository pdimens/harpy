import logging
from rich.logging import RichHandler
from rich.progress import (
    Progress,
    BarColumn,
    TaskProgressColumn,
    TextColumn,
    TimeElapsedColumn,
    SpinnerColumn,
)
from rich.console import Console
from rich.live import Live
from rich.panel import Panel
from snakemake_interface_logger_plugins.settings import OutputSettingsLoggerInterface
from snakemake_logger_plugin_harpyrich.event_handler import EventHandler


class RichLogHandler(RichHandler):
    """
    A Snakemake logger that displays job information and
    shows progress bars for rules.
    """

    def __init__(
        self,
        settings: OutputSettingsLoggerInterface,
        *args,
        **kwargs,
    ):
        self.settings = settings
        self.console = Console(log_path=False, stderr=True)
        self.progress = Progress(
            SpinnerColumn(spinner_name = "dots12", style = "blue dim", finished_text="[dim green]✓"),
            TextColumn("{task.fields[active]}", style="default"),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(bar_width=None, complete_style="yellow", finished_style="dim blue"),
            TaskProgressColumn("[progress.remaining]{task.completed}/{task.total}"),
            TimeElapsedColumn(),
            transient = True,
            auto_refresh = True,
            disable = True,
            console= self.console,
            expand=True
        )

        self.live_display = Live(
            Panel(
                self.progress,
                title="Progress",
                border_style="dim",
                padding=(0, 1, 0, 1),
            ),
            refresh_per_second=8,
            transient=True,
            console=self.console,
        )

        self.event_handler = EventHandler(
            console=self.console,
            progress=self.progress,
            live_display=self.live_display,
            dryrun=self.settings.dryrun,
            printshellcmds=self.settings.printshellcmds,
            show_failed_logs=settings.show_failed_logs,
            verbose=False
            #verbose=self.settings.verbose,
        )

        kwargs["console"] = self.console
        kwargs["show_time"] = True
        kwargs["omit_repeated_times"] = False
        kwargs["rich_tracebacks"] = True
        kwargs["tracebacks_width"] = 100
        kwargs["tracebacks_show_locals"] = False
        super().__init__(*args, **kwargs)

    def emit(self, record):
        """Process log records and delegate to event handler."""
        try:
            self.event_handler.handle(record)

        except Exception as e:
            self.handleError(
                logging.LogRecord(
                    name="RichLogHandler",
                    level=logging.ERROR,
                    pathname="",
                    lineno=0,
                    msg=f"Error in logging handler: {str(e)}",
                    args=(),
                    exc_info=None,
                )
            )

    def close(self):
        """Clean up resources."""
        self.event_handler.close()
        super().close()
