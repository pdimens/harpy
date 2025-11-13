from logging import LogRecord
from typing import Optional
from uuid import UUID
from rich.console import Console
from rich.syntax import Syntax
from rich.progress import Progress, TaskID
from rich.live import Live
from rich import box
from rich.markdown import Markdown
from rich.status import Status
from rich.panel import Panel
from rich.table import Table
import time as _time
from typing import Dict
from pathlib import Path
from snakemake_interface_logger_plugins.common import LogEvent
import snakemake_logger_plugin_harpyrich.events as events
import re
import logging


def get_time():
    return _time.strftime("%d %b %Y @ %H:%M")


def formatted_table(cols: int, left_col_style: str):
    """Convenience function that returns a table with standardized formatting"""
    _table = Table(
        show_header=False,
        pad_edge=False,
        show_edge=False,
        padding=(0, 0),
        box=box.SIMPLE,
    )
    # the column name is irrelevant b/c headers won't be shown
    _table.add_column("detail", justify="left", style=left_col_style, no_wrap=True)
    for i in range(cols - 1):
        _table.add_column(f"col_{i}", justify="left")
    return _table


def prettyprint_rule(rule: str) -> str:
    """Format the rule name to replace underscores with spaces and strip extra spaces"""
    return re.sub(" +", " ", rule.replace("_", " ")).strip()


def format_wildcards(wildcards) -> str | None:
    """Format wildcards into a string representation, if present, otherwise return None."""
    if not wildcards:
        return None

    wc_text = []
    for wc, v in wildcards.items():
        wc_text.append(f"{wc}=[cyan]{v}[/]")
    return ", ".join(wc_text)


class ProgressDisplay:
    def __init__(self, progress: Progress, live: Live, console: Console):
        self.progress = progress
        self.live_display = live
        self.rule_tasks: Dict[str, TaskID] = {}

    def add_or_update(
        self,
        rule: str,
        completed: int,
        total: int,
        visible: bool = True,
        decrement_active: bool = False,
    ):
        """
        Add a rule to the progressbar if it's not already on there, or update the progress of a rule if it is
        (which implies the job finished). Reduces the "active" field by 1 when decrement_active = True.
        """
        _rule = prettyprint_rule(rule)
        if rule not in self.rule_tasks:
            task_id = self.progress.add_task(
                description=_rule, total=total, visible=visible, active=1
            )
            self.rule_tasks[rule] = task_id
            currently_active = 1
        else:
            task_id = self.rule_tasks[rule]
            modifier = -1 if decrement_active else 1
            currently_active = self.progress.tasks[task_id].fields["active"] + modifier

        self.progress.update(
            task_id, completed=completed, total=total, active=currently_active
        )

        if completed >= total:
            self.progress.update(
                task_id, description=f"[dim]{_rule}[/]", active="[dim]-[/]"
            )

        return task_id

    def update_active(self, rule: str):
        """Increment the "active" field in a progress bar by 1"""
        task_id = self.rule_tasks[rule]
        current_task = self.progress.tasks[task_id]
        self.progress.update(task_id, active=current_task.fields["active"] + 1)

    def mark_rule_failed(self, rule: str):
        """Update progress bar for a failed rule."""
        _rule = prettyprint_rule(rule)
        if rule in self.rule_tasks:
            task_id = self.rule_tasks[rule]
            self.progress.update(
                task_id, description=f"[red]✗[/] {_rule} [red](failed)[/]"
            )

    def set_visible(self, rule: str, visible: bool = True):
        """Set visibility of a progress bar."""
        if rule in self.rule_tasks:
            task_id = self.rule_tasks[rule]
            self.progress.update(task_id, visible=visible)

    def has_tasks(self) -> bool:
        """Check if there are any active tasks."""
        return len(self.rule_tasks) > 0


class EventHandler:
    """Base class for processing Snakemake log events."""

    def __init__(
        self,
        console: Console,
        progress: Progress,
        live_display: Live,
        dryrun: bool = False,
        printshellcmds: bool = False,
        show_failed_logs: bool = False,
        verbose: bool = False,
    ):
        self.current_workflow_id: Optional[UUID] = None
        self.dryrun: bool = dryrun
        self.verbose: bool = verbose
        self.printshellcmds: bool = printshellcmds
        self.show_failed_logs: bool = show_failed_logs
        self.console = console
        self.progress = progress
        self.progress_display = ProgressDisplay(progress, live_display, self.console)
        self.jobs_info: dict = {}
        self.rule_counts: dict = {}  # {rule_name: {"total": n, "completed": m}}
        self.total_jobs = 0
        self.completed = 0
        self.conda_statuses: dict = {}  # {env_path: Status object}

    def handle(self, record: LogRecord, **kwargs) -> None:
        """Process a log record, routing to appropriate handler based on event type."""
        event_type = getattr(record, "event", None)

        if event_type:
            handler_map = {
                LogEvent.ERROR: (events.Error, self.handle_error),
                LogEvent.WORKFLOW_STARTED: (
                    events.WorkflowStarted,
                    self.handle_workflow_started,
                ),
                LogEvent.JOB_INFO: (events.JobInfo, self.handle_job_info),
                LogEvent.JOB_STARTED: (events.JobStarted, self.handle_job_started),
                LogEvent.JOB_FINISHED: (events.JobFinished, self.handle_job_finished),
                LogEvent.JOB_ERROR: (events.JobError, self.handle_job_error),
                LogEvent.SHELLCMD: (events.ShellCmd, self.handle_shellcmd),
                LogEvent.RULEGRAPH: (events.RuleGraph, self.handle_rule_graph),
                LogEvent.GROUP_INFO: (events.GroupInfo, self.handle_group_info),
                LogEvent.GROUP_ERROR: (events.GroupError, self.handle_group_error),
                LogEvent.RESOURCES_INFO: (
                    events.ResourcesInfo,
                    self.handle_resources_info,
                ),
                LogEvent.DEBUG_DAG: (events.DebugDag, self.handle_debug_dag),
                LogEvent.PROGRESS: (events.Progress, self.handle_progress),
                LogEvent.RUN_INFO: (events.RunInfo, self.handle_run_info),
            }

            handler_info = handler_map.get(event_type)
            if handler_info:
                event_class, handler_method = handler_info

                handler_method(event_class.from_record(record), **kwargs)  # type: ignore
            else:
                self.handle_generic_event(event_type, record, **kwargs)
        else:
            self.handle_generic_record(record, **kwargs)

    def handle_error(self, event_data: events.Error, **kwargs) -> None:
        """Handle error event."""
        pass

    def handle_workflow_started(
        self, event_data: events.WorkflowStarted, **kwargs
    ) -> None:
        """Handle workflow started event."""
        pass
        #self.console.rule(f"Workflow start {get_time()}", style="green")

    def handle_job_info(self, event_data: events.JobInfo, **kwargs) -> None:
        """Handle job info event with rich formatting."""
        self.jobs_info[event_data.jobid] = {
            "rule": event_data.rule_name,
            "wildcards": event_data.wildcards,
            "log": event_data.log,
        }

        self.progress_display.set_visible(event_data.rule_name, True)
        self.progress_display.update_active(event_data.rule_name)
        self.progress_display.update_active("Total Progress")

        submission_text = []
        submission_text.append(
            f"[bold light_steel_blue]◯ Started[/] {event_data.rule_name} [dim](id: {event_data.jobid})[/] [dim light_steel_blue]{get_time()}[/]"
        )
        if event_data.rule_msg:
            submission_text.append(f"[italic]{event_data.rule_msg}[/]")
        wc = format_wildcards(event_data.wildcards)
        if wc:
            submission_text.append(f"[light_steel_blue]Wildcards:[/] {wc}")
        if self.printshellcmds and self.verbose and event_data.shellcmd:
            format_cmd = re.sub(r" +", " ", event_data.shellcmd).rstrip()
            format_cmd = re.sub("^\n", "", format_cmd)
            submission_text.append("[light_steel_blue]Shell Command:[/]")
            shell_table = formatted_table(2, "default")
            cmd = Syntax(
                format_cmd,
                dedent=True,
                lexer="bash",
                tab_size=2,
                word_wrap=True,
                padding=1,
            )
            shell_table.add_row("     ", cmd)
            out_text = "\n    ".join(submission_text)
            self.console.print(out_text, shell_table, "", sep="\n")
        elif self.verbose:
            out_text = "\n    ".join(submission_text)
            self.console.print(out_text, sep="\n")

    def handle_job_started(self, event_data: events.JobStarted, **kwargs) -> None:
        """Handle job started event."""
        return

    def handle_job_finished(self, event_data: events.JobFinished, **kwargs) -> None:
        """Handle job finished event with rich formatting."""
        job_id = event_data.job_id

        if job_id in self.jobs_info:
            info = self.jobs_info[job_id]
            rule_name = info["rule"]
            wc = format_wildcards(info.get("wildcards", None))

            if rule_name in self.rule_counts:
                self.rule_counts[rule_name]["completed"] += 1
                completed = self.rule_counts[rule_name]["completed"]
                total = self.rule_counts[rule_name]["total"]

                self.progress_display.add_or_update(
                    rule_name, completed, total, decrement_active=True
                )

            self.completed += 1
            self.progress_display.add_or_update(
                "Total Progress", self.completed, self.total_jobs, decrement_active=True
            )
            if self.verbose:
                out_text = [
                    "[bold green]◉ Finished[/] "
                    + rule_name
                    + f" [dim](id: {job_id})[/] [dim green]{get_time()}[/]"
                ]
                if wc:
                    out_text.append(f"[bold green]Wildcards:[/] {wc}")
                self.console.print("\n    ".join(out_text), end="\n\n")

    def handle_shellcmd(self, event_data: events.ShellCmd, **kwargs) -> None:
        """Handle shell command event with syntax highlighting."""
        return

    def handle_job_error(self, event_data: events.JobError, **kwargs) -> None:
        """Handle job error event."""
        job_id = event_data.jobid
        if job_id in self.jobs_info:
            info = self.jobs_info[job_id]
            rule_name = info["rule"]
            wc = format_wildcards(info["wildcards"])
            failed_text = (
                f"[bold yellow]✗ Failed[/] {rule_name} [dim yellow](id: {job_id})[/]"
            )
            if wc:
                failed_text += f"\n    [bold yellow]Wildcards:[/] {wc}"
            self.console.log(failed_text)
            # TODO Not working like it's supposed to
            # self.console.log(info["log"])
            # if self.show_failed_logs:
            for _log in info["log"]:
                self.console.rule(f"[bold]Log file: {_log}", style="yellow")
                self.console.print(Path(_log).read_text(), highlight=False)
        else:
            self.console.log(f"[bold yellow]✗ Failed job_id: {job_id})[/]")

    def handle_group_info(self, event_data: events.GroupInfo, **kwargs) -> None:
        """Handle group info event."""
        pass

    def handle_group_error(self, event_data: events.GroupError, **kwargs) -> None:
        """Handle group error event."""
        pass

    def handle_resources_info(self, event_data: events.ResourcesInfo, **kwargs) -> None:
        """Handle resources info event."""
        pass

    def handle_debug_dag(self, event_data: events.DebugDag, **kwargs) -> None:
        """Handle debug DAG event."""
        pass

    def handle_progress(self, event_data: events.Progress, **kwargs) -> None:
        """Handle progress event."""
        pass

    def handle_rule_graph(self, event_data: events.RuleGraph, **kwargs) -> None:
        """Handle rule graph event."""
        pass

    def handle_run_info(self, event_data: events.RunInfo, **kwargs) -> None:
        """Handle run info event - sets up progress bars."""
        try:
            self.dag_status.stop()
        except NameError:
            pass
        self.total_jobs = event_data.total_job_count

        if self.total_jobs > 0:
            self.total_progress_task = self.progress_display.add_or_update(
                "Total Progress", 0, self.total_jobs
            )
            # self.console.log(f"Processing Workflow: {self.total_jobs} jobs", style="blue")
            self.progress.disable = False
            # end any existing conda statuses
            for status in self.conda_statuses.values():
                status.stop()
            self.conda_statuses.clear()
            self.progress_display.live_display.start()

        for rule, count in event_data.per_rule_job_counts.items():
            if count > 0:
                self.rule_counts[rule] = {"total": count, "completed": 0}
                self.progress_display.add_or_update(rule, 0, count, visible=False)

    def handle_generic_event(
        self, event_type: LogEvent, record: LogRecord, **kwargs
    ) -> None:
        """Handle events that don't have a specific handler defined."""
        pass

    def handle_generic_record(self, record: LogRecord, **kwargs) -> None:
        """Handle log records that don't have an event type."""
        message = record.getMessage()

        if not self.should_log_message(record, message):
            return
        conda_depwarn = "Your conda installation is not configured to use" in message
        if conda_depwarn and self.verbose:
            self.console.print(
                Panel(
                    Markdown(
                        "Adding `defaults` to the conda channel list implicitly is deprecated. To fix this, read [this guide](https://conda-forge.org/docs/user/tipsandtricks.html)."
                    ),
                    title="Warning: conda channel configuration",
                    border_style="yellow",
                )
            )
            return

        # Check for conda environment creation start
        conda_create_match = re.search(
            r"Creating conda environment (.+?)\.\.\..*", message
        )
        if conda_create_match:
            env_path = conda_create_match.group(1)
            env_name = Path(env_path).name
            self._start_conda_status(env_name)
            return

        # Check for conda environment creation completion
        conda_done_match = re.search(
            r"Environment for (.+?) created \(location: (.+?)\)", message
        )
        if conda_done_match:
            env_path = conda_done_match.group(1)
            env_name = Path(env_path).name
            self._complete_conda_status(env_name)
            return

        if "Complete log" in message:
            return

        build_dag_match = re.search("^Building DAG of jobs", message)
        if build_dag_match:
            self._start_dag_status()

    def _start_dag_status(self):
        """Start a spinning status for conda environment creation."""
        status = Status(
            "Building workflow graph...",
            console=self.console,
            spinner="dots",
        )
        status.start()
        self.dag_status = status

    def _start_conda_status(self, env_name: str):
        """Start a spinning status for conda environment creation."""
        try:
            self.dag_status.stop()
        except NameError:
            pass
        status = Status(
            f"Creating conda environment [cyan]{env_name}[/cyan]...",
            console=self.console,
            spinner="dots",
        )
        status.start()
        self.conda_statuses[env_name] = status

    def _complete_conda_status(self, env_name: str):
        """Complete the conda environment creation status."""
        if env_name in self.conda_statuses:
            status = self.conda_statuses[env_name]
            if self.verbose:
                self.console.log(
                    f"[green]◉ Created[/] conda environment [cyan]{env_name}[/cyan]"
                )
            status.stop()
            del self.conda_statuses[env_name]

        else:
            if self.verbose:
                self.console.log(
                    f"[green]◉ Created[/] conda environment [cyan]{env_name}[/cyan]"
                )
            return

    def close(self):
        """Clean up any active statuses."""
        for status in self.conda_statuses.values():
            status.stop()
        self.conda_statuses.clear()
        self.progress_display.live_display.stop()
        self.progress_display.progress.stop()

    def should_log_message(self, record, message):
        """Determine if we should log this message based on content filtering."""

        if message == "None":
            return False

        if record.levelno >= logging.ERROR:
            return True

        skip_patterns = [
            "^Select jobs to execute",
            "^Assuming unrestricted shared filesystem",
            "^Using shell:",
            "^host:",
            "^Provided cores:",
            "^Rules claiming more threads will be",
            r"^Execute \d+ jobs.",
            # "^Building DAG of jobs.",
            "^Activating conda env",
        ]

        for pattern in skip_patterns:
            if re.search(pattern, message):
                return False

        return True
