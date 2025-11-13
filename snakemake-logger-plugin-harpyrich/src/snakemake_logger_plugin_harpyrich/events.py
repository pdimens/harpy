import uuid
from dataclasses import dataclass, field
from logging import LogRecord
from typing import Any, Dict, List, Optional

# redefine these here until they are in logger interface


@dataclass
class Error:
    exception: Optional[str] = None
    location: Optional[str] = None
    rule: Optional[str] = None
    traceback: Optional[str] = None
    file: Optional[str] = None
    line: Optional[str] = None

    @classmethod
    def from_record(cls, record: LogRecord) -> "Error":
        return cls(
            exception=getattr(record, "exception", None),
            location=getattr(record, "location", None),
            rule=getattr(record, "rule", None),
            traceback=getattr(record, "traceback", None),
            file=getattr(record, "file", None),
            line=getattr(record, "line", None),
        )


@dataclass
class WorkflowStarted:
    workflow_id: uuid.UUID
    snakefile: Optional[str]

    def __post_init__(self) -> None:
        if not isinstance(self.snakefile, str):
            try:
                # Try to convert to string - this should work for PosixPath and other path-like objects
                self.snakefile = str(self.snakefile)
            except (TypeError, ValueError) as e:
                raise ValueError(f"Could not convert snakefile to string: {e}")

    @classmethod
    def from_record(cls, record: LogRecord) -> "WorkflowStarted":
        return cls(
            workflow_id=getattr(record, "workflow_id"),
            snakefile=getattr(record, "snakefile", None),
        )


@dataclass
class JobInfo:
    jobid: int
    rule_name: str
    threads: int
    input: List[str] = field(default_factory=list)
    output: List[str] = field(default_factory=list)
    log: List[str] = field(default_factory=list)
    benchmark: List[str] = field(default_factory=list)
    rule_msg: Optional[str] = None
    wildcards: Dict[str, Any] = field(default_factory=dict)
    reason: Optional[str] = None
    shellcmd: Optional[str] = None
    priority: Optional[int] = None
    resources: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_record(cls, record: LogRecord) -> "JobInfo":
        resources = {}
        if hasattr(record, "resources") and hasattr(record.resources, "_names"):  # type: ignore
            resources = {
                name: value
                for name, value in zip(record.resources._names, record.resources)  # type: ignore
                if name not in {"_cores", "_nodes"}
            }

        return cls(
            jobid=getattr(record, "jobid", 0),
            rule_name=getattr(record, "rule_name", ""),
            threads=getattr(record, "threads", 1),
            rule_msg=getattr(record, "rule_msg", None),
            wildcards=getattr(record, "wildcards", {}),
            reason=getattr(record, "reason", None),
            shellcmd=getattr(record, "shellcmd", None),
            priority=getattr(record, "priority", None),
            input=getattr(record, "input", []),
            log=getattr(record, "log", []),
            output=getattr(record, "output", []),
            benchmark=getattr(record, "benchmark", []),
            resources=resources,
        )


@dataclass
class JobStarted:
    job_ids: List[int] = field(default_factory=list)

    @classmethod
    def from_record(cls, record: LogRecord) -> "JobStarted":
        jobs = getattr(record, "jobs", [])

        if jobs is None:
            jobs = []
        elif isinstance(jobs, int):
            jobs = [jobs]

        return cls(job_ids=jobs)


@dataclass
class JobFinished:
    job_id: int

    @classmethod
    def from_record(cls, record: LogRecord) -> "JobFinished":
        return cls(job_id=getattr(record, "job_id"))


@dataclass
class ShellCmd:
    jobid: int
    shellcmd: Optional[str] = None
    rule_name: Optional[str] = None

    @classmethod
    def from_record(cls, record: LogRecord) -> "ShellCmd":
        return cls(
            jobid=getattr(record, "jobid", 0),
            shellcmd=getattr(record, "shellcmd", ""),
            rule_name=getattr(record, "name", None),
        )


@dataclass
class JobError:
    jobid: int

    @classmethod
    def from_record(cls, record: LogRecord) -> "JobError":
        return cls(
            jobid=getattr(record, "jobid", 0),
        )


@dataclass
class GroupInfo:
    group_id: int
    jobs: List[Any] = field(default_factory=list)

    @classmethod
    def from_record(cls, record: LogRecord) -> "GroupInfo":
        return cls(
            group_id=getattr(record, "group_id", 0), jobs=getattr(record, "jobs", [])
        )


@dataclass
class GroupError:
    groupid: int
    aux_logs: List[Any] = field(default_factory=list)
    job_error_info: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_record(cls, record: LogRecord) -> "GroupError":
        return cls(
            groupid=getattr(record, "groupid", 0),
            aux_logs=getattr(record, "aux_logs", []),
            job_error_info=getattr(record, "job_error_info", {}),
        )


@dataclass
class ResourcesInfo:
    nodes: List[str] = field(default_factory=list)
    cores: Optional[int] = None
    provided_resources: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_record(cls, record: LogRecord) -> "ResourcesInfo":
        if hasattr(record, "nodes"):
            return cls(nodes=getattr(record, "nodes", []))
        elif hasattr(record, "cores"):
            return cls(cores=record.cores)  # type: ignore
        elif hasattr(record, "provided_resources"):
            return cls(provided_resources=getattr(record, "provided_resources", {}))
        else:
            return cls()


@dataclass
class DebugDag:
    status: Optional[str] = None
    job: Optional[Any] = None
    file: Optional[str] = None
    exception: Optional[str] = None

    @classmethod
    def from_record(cls, record: LogRecord) -> "DebugDag":
        return cls(
            status=getattr(record, "status", None),
            job=getattr(record, "job", None),
            file=getattr(record, "file", None),
            exception=getattr(record, "exception", None),
        )


@dataclass
class Progress:
    done: int
    total: int

    @classmethod
    def from_record(cls, record: LogRecord) -> "Progress":
        return cls(done=getattr(record, "done", 0), total=getattr(record, "total", 0))


@dataclass
class RuleGraph:
    rulegraph: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_record(cls, record: LogRecord) -> "RuleGraph":
        return cls(rulegraph=getattr(record, "rulegraph", {}))


@dataclass
class RunInfo:
    per_rule_job_counts: Dict[str, int] = field(default_factory=dict)
    total_job_count: int = 0

    @classmethod
    def from_record(cls, record: LogRecord) -> "RunInfo":
        all_stats = getattr(record, "stats", {})

        per_rule_job_counts = {k: v for k, v in all_stats.items() if k != "total"}

        total_job_count = all_stats.get("total", 0)
        return cls(
            per_rule_job_counts=per_rule_job_counts, total_job_count=total_job_count
        )
