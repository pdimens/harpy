
from harpy.common.printing import HarpyPrint
from dataclasses import dataclass, field
import re
import sys

from rich.markup import escape


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
    logs: list[str] = field(default_factory=list)
    log_contents: dict[str, str] = field(default_factory=dict)
    resources: list[str] = field(default_factory=list)

_HEADER_RE = re.compile(r"^(?:Error in rule|Error in group|rule) (\w+):\s*$")
_KEY_RE = re.compile(r"^\s{4}(\S[\w -]*):\s?(.*)$")
_LOGFILE_HEADER_RE = re.compile(r"^Logfile (\S+)(?: \(send to storage\))?:\s*$")
_LOGFILE_NOTFOUND_RE = re.compile(r"^Logfile \S+.*not found\.\s*$")
_ERROR_RULE_RE = re.compile(
    r"^\s*Error in rule (\w+):\s*$"
)
_KEY_NAMES = {
    "jobid", "input", "output", "log", "message", "conda-env",
    "shell", "external_jobid", "resources", "threads", "jobs",
}
_KEY_RE = re.compile(
    r"^\s*(" + "|".join(re.escape(k) for k in _KEY_NAMES) + r"):\s?(.*)$"
)
_HEADER_RE = re.compile(r"^\s*(?:Error in rule|Error in group|rule) (\w+):\s*$")
_GROUP_MSG_RE = re.compile(r"^Error in group (\S+)$")
_LATENCY_TRIGGER_RE = re.compile(r"--latency-wait:\s*$")
_CORRUPTED_TRIGGER_RE = re.compile(r"corrupted:\s*$")


def is_log_sep(line: str) -> bool:
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


def fmt_missing_block(text: str) -> str:
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

def fmt_shell_cmd(text: str, indent: str = "    ") -> str:
    """
    Reindent a shell block by brace nesting, not by whatever whitespace
    snakemake printed (which drops/adds blank lines and dedents inconsistently).
    """
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



class ErrorHandler():
    def __init__(self, err):
        self.errortext = _Pushback(iter(err))
        self.rules: list[SnakeRule] = []
        self.missingoutput: list[str] = []
        # printing config
        self.hp= HarpyPrint()
        self.hp.console.tab_size = 4
        self.hp.console._highlight = False
        self.hp.console.soft_wrap = True

    def process(self):
        '''
        final processing of the snakemake stderr text after an error has occured,
        returns early if ongoing or successful exit, otherwise processess the error text.
        '''
        # ---------- shortcut to FileNotFoundError
        line = next(self.errortext)
        if line.strip().startswith("FileNotFound"):
            if "envs/" in line and ".yaml'" in line:
                self.hp.print("[red]Missing conda environment yaml file:[/][yellow]\n  " + line.split(":")[-1].replace("'", ""))
            else:
                self.hp.print(line.strip(), style = "red")
            return

        # ---------- pick out conda-env errors
        if "Could not create conda" in line:
            for i in self.errortext:
                if "To search for alternate" in i:
                    break
                if i.strip():
                    if i.lstrip().startswith("-"):
                        self.hp.print(i.rstrip(), soft_wrap = True, width = 2000, style = "bold red")
                    else:
                        self.hp.print(i.rstrip(), soft_wrap = True, width = 2000, style = "red")
            return

        if ("Error" in line or "Exception" in line or "Missing input files" in line) and not ("RuleException" in line or "CalledProcessError" in line):
            #self.rule("[bold]Source of Error", style = "dim")
            self.hp.print(line, highlight=False, soft_wrap = True, end = "", style = "red")
            for i in self.errortext:
                self.hp.print(i, highlight = False, soft_wrap = True, end = "", style="red")
            return

        if "but some output files are missing" in line:
            self.missingoutput.append(line)
            for i in self.errortext:
                if "Shutting down, this might" in i:
                    break
                elif "but some output files are missing" not in i:
                    self.missingoutput[-1] += i
                else:
                    self.missingoutput.append(i)

        for i in self.errortext:
            if "Exiting because a job execution failed. Look below for error messages" in i:
                break

        for i in self.errortext:
            if not line.strip():
                continue

            if _ERROR_RULE_RE.match(line):
                self.errortext.push(line)
                break

            if line.lstrip().startswith("Logfile "):
                self.errortext.push(line)
                break
            if "(100%) done" in i:
                break

            if "RuleException" in i:
                sys.exit(1)

            m = _KEY_RE.match(i)
            if m and m.group(1) == "message":
                gm = _GROUP_MSG_RE.match(m.group(2).strip())
                if gm:
                    self._skip_group_block()
                    continue

            hm = _ERROR_RULE_RE.match(i.strip())
            if hm:
                rule = self._parse_rule_block()
                rule.name = hm.group(1)
                self.rules.append(rule)
                print(rule)

            elif i.startswith("Complete log"):
                break

            elif i.startswith("WorkflowError"):
                break

        if self.missingoutput:
            self.hp.print("\n[bold dim]──── ⚠ Error Reported by Snakemake")
            self.hp.print(fmt_missing_block(self.missingoutput[0]), style = "red")
        elif len(self.rules) > 1:
            grp = " → ".join(rule.name for rule in self.rules)
            self.hp.rule(
                f"[default bold]Triggering Group[/][yellow bold] {grp}",
                style="yellow"
            )
            for rule in self.rules:
                self.hp.print(f"\n──── {rule.name}", style = "bold yellow")
                self.print(rule)

        elif len(self.rules) == 1:
            self.hp.rule(
                f"[default bold]Triggering Rule[/][yellow bold] {self.rules[0].name}",
                style="yellow"
            )
            self.print(self.rules[0])    

    def _parse_rule_block(self) -> SnakeRule:
        """Parse one complete `Error in rule` block."""
        rule = SnakeRule()
        key = None

        for line in self.errortext:
            if not line.strip():
                continue

            # A new rule / end-of-rule marker.
            if _ERROR_RULE_RE.match(line):
                self.errortext.push(line)
                break

            if line.lstrip().startswith("Logfile "):
                self.errortext.push(line)
                break

            if line.lstrip().startswith("Complete log"):
                self.errortext.push(line)
                break

            if line.lstrip().startswith("WorkflowError"):
                self.errortext.push(line)
                break

            m = _KEY_RE.match(line)
            if m:
                key, val = m.group(1), m.group(2)
                if key == "jobid":
                    rule.jobid = int(val)
                elif key == "input":
                    rule.input = val
                elif key == "output":
                    rule.output = val
                elif key == "log":
                    pass
                elif key == "message":
                    rule.message = val
                elif key == "conda-env":
                    rule.env = val
                elif key == "resources":
                    rule.resources = [r.strip() for r in val.split(",")]
                elif key == "shell":
                    rule.cmd = ""
            elif key == "shell" and line[:1].isspace():
                rule.cmd += line + "\n"

            elif line[:1].isspace():
                continue

            else:
                self.errortext.push(line)
                break

        self._consume_logfile_blocks(rule)
        rule.cmd = fmt_shell_cmd(rule.cmd)

        return rule

    def _consume_logfile_blocks(self, rule: SnakeRule):
        """Snakemake prints 'Logfile <path>:' + fenced content, or
        'Logfile <path> ... not found.', right after a rule's error block.
        Consumes zero or more of these, pushing back the first line that isn't one."""
        for line in self.errortext:
            m = _LOGFILE_HEADER_RE.match(line.strip())
            if m:
                content = []
                for j in self.errortext:
                    if is_log_sep(j):
                        if content:
                            break  # closing fence
                        continue  # opening fence
                    content.append(j)

                logtext = escape(re.sub(r'\n{3,}', '\n\n', "".join(content)).removeprefix("    "))
                # purge out all unnecessary papermill error text
                if "papermill.exceptions.PapermillExecutionError:" in logtext:
                    logtext = logtext.partition("papermill.exceptions.PapermillExecutionError:")[-1]
                    chunks = logtext.split("\n\n")
                    filtered = [c for c in chunks if not c.startswith("File ")]
                    logtext = "\n\n".join(filtered)
                rule.logs.append(m.group(1))
                rule.log_contents[m.group(1)] = logtext.strip()
                continue
            if _LOGFILE_NOTFOUND_RE.match(line.strip()):
                continue  # path already in rule.logs, nothing to capture
            self.errortext.push(line)
            break

    def _parse_rule_stub(self) -> SnakeRule:
        """Manifest entry inside a group's 'jobs:' list — jobid/output/log only,
        no shell/input/env. Ends at next 'rule X:' header or first unrecognized line."""
        rule = SnakeRule()
        for line in self.errortext:
            if not line.strip():
                continue
            if _HEADER_RE.match(line):
                self.errortext.push(line)
                break
            m = _KEY_RE.match(line)
            if not m or m.group(1) not in ("jobid", "output", "log"):
                self.errortext.push(line)
                break
            key, val = m.group(1), m.group(2)
            if key == "jobid":
                rule.jobid = int(val)
            elif key == "output":
                rule.output = val
            elif key == "log":
                rule.logfile = val
                rule.logs.append(val)
        return rule


    def _skip_group_block(self):
        """Skip the incomplete group jobs manifest."""
        for line in self.errortext:
            if _ERROR_RULE_RE.match(line):
                self.errortext.push(line)
                return

#    def _parse_group_block(self, groupid: str) -> SnakeGroup:
#        """
#        Called right after a 'message: Error in group <id>' line. Consumes the
#        'jobs:' manifest of rule stubs that follows.
#        """
#        group = SnakeGroup(groupid=groupid)
#        for line in self.errortext:
#            if not line.strip():
#                continue
#            m = _KEY_RE.match(line)
#            if m and m.group(1) == "jobs":
#                continue  # bare 'jobs:' line, entries follow
#            hm = _HEADER_RE.match(line)
#            if hm:
#                stub = self._parse_rule_stub()
#                stub.name = hm.group(1)
#                group.rules.append(stub)
#                continue
#            self.errortext.push(line)
#            break
#        return group
#
    def print(self, rule: SnakeRule):
        'Print a nicely formatted Snakemake rule error'
        self.hp.print("input:", style = "bold default")
        self.hp.print("  " + rule.input.replace(", ", "\n  "), style = "red")

        self.hp.print("output:", style = "bold default")
        self.hp.print("  " + rule.output.replace(", ", "\n  "), style = "red")

        if rule.logs:
            self.hp.print("log:", style = "bold default")
            for i in rule.logs:
                _i = i.replace("(check log file(s) for error details)", "").strip()
                self.hp.print("  " + _i, style = "red")

        if rule.env:
            env_clean = rule.env.split("/")[:-2]
            self.hp.print(
                f"[bold default]conda-env:[/] [red]{'/'.join(env_clean)}[/]"
            )

        if rule.message != "None":
            self.hp.print("Message:", style = "bold default")
            self.hp.print(rule.message, style = "red")

        if rule.cmd:
            self.hp.print("")
            self.hp.print("[bold dim]──── ❯ Command Invoked")
            self.hp.shell(rule.cmd)

        if rule.logs:
            self.hp.print("")

        #for name, content in rule.logs.items():
        for i in rule.logs:
            self.hp.print(f"──── 🗎 {i}", style = "bold dim")
            self.hp.print(rule.log_contents[i], style = "red")

