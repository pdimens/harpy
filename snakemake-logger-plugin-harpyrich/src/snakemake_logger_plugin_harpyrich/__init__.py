from snakemake_interface_logger_plugins.base import LogHandlerBase
from snakemake_logger_plugin_harpyrich.handler import (
    RichLogHandler,
)
from rich.console import Console


class LogHandler(LogHandlerBase, RichLogHandler):  # type: ignore
    def __post_init__(self) -> None:
        """
        Any additional setup after initialization.
        """
        console = Console(
            stderr=not self.common_settings.stdout,
        )
        RichLogHandler.__init__(self, console=console, settings=self.common_settings)

    def emit(self, record):
        """Delegate emit to RichLogHandler"""
        return RichLogHandler.emit(self, record)

    @property
    def writes_to_stream(self) -> bool:
        """
        Whether this plugin writes to stderr/stdout
        """
        return True

    @property
    def writes_to_file(self) -> bool:
        """
        Whether this plugin writes to a file
        """
        return False

    @property
    def has_filter(self) -> bool:
        """
        Whether this plugin attaches its own filter
        """
        return True

    @property
    def has_formatter(self) -> bool:
        """
        Whether this plugin attaches its own formatter
        """
        return True

    @property
    def needs_rulegraph(self) -> bool:
        """
        Whether this plugin requires the DAG rulegraph.
        """
        return False
