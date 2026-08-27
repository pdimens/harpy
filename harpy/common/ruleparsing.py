
import re
from dataclasses import dataclass, field

@dataclass
class SnakeRule:
    '''
    Dataclass to hold snakemake rule information for error printing
    '''
    name: str = ""
    jobid: int = 0
    input: str = ""
    output: str = ""
    message: str = ""
    env: str = ""
    cmd: str = ""
    logs: dict[str, str] = field(default_factory=dict)
    resources: list[str] = field(default_factory=list)

_HEADER_RE = re.compile(r"^(?:Error in rule|Error in group|rule) (\w+):\s*$")
_KEY_RE = re.compile(r"^\s{4}(\S[\w -]*):\s?(.*)$")
_LOGFILE_HEADER_RE = re.compile(r"^Logfile (\S+)(?: \(send to storage\))?:\s*$")
_LOGFILE_NOTFOUND_RE = re.compile(r"^Logfile \S+.*not found\.\s*$")

def _is_separator(line: str) -> bool:
    '''
    Identify the ======== lines that denote a log file's contents in a snakemake rule error
    '''
    s = line.strip()
    return len(s) >= 3 and set(s) == {"="}

class _Pushback:
    """
    Wraps a single-pass line iterator with one-line lookahead/pushback.
    This is used for processing snakemake output for error handling.
    """
    __slots__ = ("_it", "_buf")
    def __init__(self, it):
        self._it = it
        self._buf = []
    def __iter__(self):
        return self
    def __next__(self):
        return self._buf.pop() if self._buf else next(self._it)
    def push(self, line):
        self._buf.append(line)

_LATENCY_TRIGGER_RE = re.compile(r"--latency-wait:\s*$")
_CORRUPTED_TRIGGER_RE = re.compile(r"corrupted:\s*$")

def format_missing_output_block(text: str) -> str:
    '''
    Parsing for the MissingOutputException text from snakemake to enable sensible indentation. Returns
    the nicely-formatted text.
    '''
    lines = text.splitlines()
    out = []
    i = 0
    while i < len(lines):
        line = lines[i]
        out.append(line)
        if _LATENCY_TRIGGER_RE.search(line):
            i += 1
            while i < len(lines) and lines[i].strip() and not _CORRUPTED_TRIGGER_RE.search(lines[i]):
                out.append("  " + lines[i].strip())
                i += 1
            continue
        if _CORRUPTED_TRIGGER_RE.search(line):
            i += 1
            if i < len(lines) and lines[i].strip():
                files = [f.strip() for f in lines[i].split(",")]
                out.extend("  " + f for f in files)
                i += 1
            continue
        i += 1
    return "\n".join(out)

def format_shell_cmd(text: str, indent: str = "    ") -> str:
    """Reindent a shell block by brace nesting, not by whatever whitespace
    snakemake printed (which drops/adds blank lines and dedents inconsistently)."""
    out = []
    depth = 0
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue  # drop snakemake's injected blank lines
        dedent_now = s.startswith("}")
        cur_depth = max(depth - 1, 0) if dedent_now else depth
        out.append(indent * cur_depth + s)
        depth = max(depth + s.count("{") - s.count("}"), 0)
    return "\n".join(out)