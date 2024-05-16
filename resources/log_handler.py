"""Handler for Snakemake log messages."""
import os
import re
import sys
import time
from pathlib import Path
from datetime import timedelta

from alive_progress import alive_bar
from snakemake import logger

from sparv.core import paths


class LogHandler:
    """Class providing a log handler for Snakemake."""

    # The second unicode character below is a zero width space, needed for the progressbar to calculate the correct
    # width. It uses len() which doesn't take into account that the bird character takes up two columns when printed.
    icon = "\U0001f426\U0000200b"

    def __init__(self, progressbar=False, summary=False):
        """Initialize log handler.

        Args:
            progressbar: Set to True to enable progress bar. Disabled by default.
            summary: Set to True to write a final summary (elapsed time). Disabled by default.
        """
        self.use_progressbar = progressbar
        self.show_summary = summary
        self.finished = False
        self.messages = []
        self.real_errors = []
        self.load_errors = None
        self.handled_load_errors = set()
        self.start_time = time.time()
        self.exit_message = None

        # Progress bar related variables
        self.bar_mgr = None
        self.exit = lambda *x: None
        self.bar = None
        self.last_percentage = 0

    def setup_bar(self, total: int):
        """Initialize the progress bar."""
        print()
        self.bar_mgr = (alive_bar(total, enrich_print=False, title=LogHandler.icon, length=30))
        self.exit = type(self.bar_mgr).__exit__
        self.bar = type(self.bar_mgr).__enter__(self.bar_mgr)

    def log_handler(self, msg):
        """Log handler for Snakemake displaying a progress bar."""
        level = msg["level"]

        if level == "run_info" and self.use_progressbar:
            if self.bar is None:
                # Get number of jobs
                jobs: str = msg["msg"].split("\t")[-1]
                if jobs.isdigit():
                    self.setup_bar(int(jobs))

        elif level == "progress":
            if self.use_progressbar:
                # Set up progress bar if needed
                if self.bar is None:
                    self.setup_bar(msg["total"])

                # Advance progress
                self.bar()

                # Print regular updates if output is not a terminal (i.e. doesn't support the progress bar)
                if not sys.stdout.isatty():
                    percentage = (100 * msg["done"]) // msg["total"]
                    if percentage > self.last_percentage:
                        self.last_percentage = percentage
                        print(f"Progress: {percentage}%")

            if msg["done"] == msg["total"]:
                self.stop()

        elif level == "job_info" and self.use_progressbar:
            if msg["msg"] and self.bar is not None:
                if msg["msg"].startswith("EXIT_MESSAGE: "):
                    self.exit_message = msg["msg"][14:]
                else:
                    # Only update status message, don't advance progress
                    self.bar.text(msg["msg"])

        elif level == "info":
            if msg["msg"] == "Nothing to be done.":
                logger.text_handler(msg)

        elif level == "error" and ("exit status 123" in msg["msg"] or (
            "SystemExit" in msg["msg"] and "123" in msg["msg"]
        )):
            # Exit status 123 means a SparvErrorMessage exception
            # Find log files tied to this process
            for log_file in paths.log_dir.glob("{}.*".format(os.getpid())):
                log_info = re.match(r"[^.]+\.[^.]+\.([^.]+)\.([^.]+)\.([^.]+)", log_file.stem)
                self.messages.append((log_info.groups(), log_file.read_text()))
                log_file.unlink()

        elif level in ("warning", "error", "job_error"):
            # Save other errors and warnings for later
            self.real_errors.append(msg)

        elif level == "dag_debug" and "job" in msg:
            # If a module can't be used due to missing config variables, a log is created. Read those log files and
            # save to a dictionary.
            if self.load_errors is None:
                self.load_errors = {}
                for log_file in paths.log_dir.glob("{}.load_error.*".format(os.getpid())):
                    annotator_name = log_file.stem.split(".", 2)[2].replace(".", "::")
                    self.load_errors.setdefault(annotator_name, [])
                    self.load_errors[annotator_name].append(log_file.read_text())
                    log_file.unlink()

            # Check the rules used by the current operation, and see if any is unusable
            if msg["status"] == "candidate":
                job_name = str(msg["job"])
                if job_name in self.load_errors and job_name not in self.handled_load_errors:
                    self.handled_load_errors.add(job_name)
                    for error_message in self.load_errors[job_name]:
                        self.messages.append((["error"] + job_name.split("::"), error_message))

    def stop(self):
        """Stop the progress bar and output any error messages."""
        # Make sure this is only run once
        if not self.finished:
            # Stop progress bar
            if self.bar is not None:
                self.exit(self.bar_mgr, None, None, None)

            self.finished = True

            # Remove Snakemake log file if empty (which it will be unless errors occurred)
            snakemake_log_file = logger.get_logfile()
            if snakemake_log_file is not None:
                log_file = Path(snakemake_log_file)
                if log_file.is_file() and log_file.stat().st_size == 0:
                    log_file.unlink()

            print()

            # Print user-friendly error messages
            if self.messages:
                logger.logger.error("Job execution failed with the following message{}:".format(
                    "s" if len(self.messages) > 1 else ""))
                for message in self.messages:
                    (_message_type, module_name, f_name), msg = message
                    logger.logger.error("\n[{}:{}]\n{}".format(module_name, f_name, msg))
            # Defer to Snakemake's default log handler for other errors
            elif self.real_errors:
                for error in self.real_errors:
                    logger.text_handler(error)

            if self.exit_message:
                logger.logger.info(self.exit_message)

            if self.show_summary:
                if self.messages or self.real_errors:
                    print()
                elapsed = round(time.time() - self.start_time)
                logger.logger.info("Time elapsed: {}".format(timedelta(seconds=elapsed)))


def exit_with_message(error, snakemake_pid, pid, module_name, function_name):
    """Save error message to temporary file (to be read by log handler) and exit."""
    log_file = paths.log_dir / "{}.{}.error.{}.{}.log".format(snakemake_pid, pid or 0, module_name, function_name)
    log_file.parent.mkdir(parents=True, exist_ok=True)
    log_file.write_text(str(error))
    sys.exit(123)