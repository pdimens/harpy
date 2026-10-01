
from harpy.common.printing import HarpyPrint
from dataclasses import dataclass, field
import os
import re
import sys

from rich.markup import escape


@dataclass
class SnakeRule:
    '''Dataclass to hold snakemake rule information for error printing'''
    name: str = ''
    jobid: int = 0
    input: str = ''
    output: str = ''
    message: str = ''
    env: str = ''
    cmd: str = ''
    apptainerprefix: str = ''
    logs: list[str] = field(default_factory=list)
    log_contents: dict[str, str] = field(default_factory=dict)
    resources: list[str] = field(default_factory=list)

#_HEADER_RE = re.compile(r"^(?:Error in rule|Error in group|rule) (\w+):\s*$")
_KEY_RE = re.compile(r"^\s{4}(\S[\w -]*):\s?(.*)$")
_ERROR_RULE_RE = re.compile(r"^\s*Error in rule (\w+):\s*$")
_HEADER_RE = re.compile(r"^\s*(?:Error in rule|Error in group|rule) (\w+):\s*$")
_GROUP_MSG_RE = re.compile(r"^Error in group (\S+)$")
_LOGFILE_HEADER_RE = re.compile(r"^Logfile (\S+)(?: \(send to storage\))?:\s*$")
_LOGFILE_NOTFOUND_RE = re.compile(r"^Logfile \S+.*not found\.\s*$")
_LOGFILE_NOTFOUND_PATH_RE = re.compile(r"^Logfile (\S+).*not found\.\s*$")
_LOGFILE_EMPTY_RE = re.compile(r"^Logfile (\S+): empty file\s*$")
_LOGFILE_BINARY_RE = re.compile(r"^Logfile (\S+) is not a text file\.\s*$")
_LOG_LIST_RE = re.compile(r"^\s{4}log:\s*(.+?)(?:\s*\(check log file\(s\) for error details\))?\s*$")
_MESSAGE_RE = re.compile(r"^\s{4}message:\s*(.*)$")
_EXTERNAL_JOBID_RE = re.compile(r"^\s{4}external_jobid:\s*(.+?)\s*$")
_LATENCY_TRIGGER_RE = re.compile(r"--latency-wait:\s*$")
_CORRUPTED_TRIGGER_RE = re.compile(r"corrupted:\s*$")
#_EXITING_RE = re.compile(r"Exiting because a job execution failed\. Look ")
_KEY_NAMES = {
    'jobid', 'input', 'output', 'log', 'message', 'conda-env',
    'shell', 'external_jobid', 'resources', 'threads', 'jobs',
}
_KEY_RE = re.compile(
    r"^\s*(" + "|".join(re.escape(k) for k in _KEY_NAMES) + r"):\s?(.*)$"
)

def is_log_sep(line: str) -> bool:
    '''
    Identify the ======== lines that denote a log file's contents in a snakemake rule error
    '''
    s = line.strip()
    return len(s) >= 3 and set(s) == {'='}

@dataclass
class HPCErrors:
    '''What could be scraped from the snakemake output of a failed workflow run with a scheduler (executor) plugin'''
    rules: list[str] = field(default_factory=list)
    messages: list[str] = field(default_factory=list)
    jobids: list[str] = field(default_factory=list)
    # log file path -> contents, or None if snakemake reported it couldn't read the file
    logs: dict[str, str | None] = field(default_factory=dict)
    # every log file path listed on a `log:` line (rule logs, then scheduler logs), in order of appearance
    listed: list[str] = field(default_factory=list)

    def __bool__(self) -> bool:
        return bool(self.rules or self.messages or self.jobids or self.logs or self.listed)


def scrape_hpc_errors(lines) -> HPCErrors:
    '''
    Best-effort extraction of the useful bits of snakemake's error output when jobs were run by a
    scheduler plugin (slurm, lsf, googlebatch, cluster-generic, ...). The output of these is too varied and
    interleaved to parse as a structured rule block, so this doesn't try: it only picks out the things that
    look the same regardless of the plugin, ignoring everything else (including where they appear):

    - `Error in rule X:` names and `message:` lines
    - `external_jobid:` (the scheduler's job ID)
    - `Logfile PATH:` blocks, where snakemake prints the contents of the failed job's log files. These
      include the scheduler's own log when the plugin provides one, which plugins do by handing it to snakemake
      as an auxiliary log.
    - `log:` lines, for the paths of every log file involved. Snakemake stops printing log contents at the first
      file it can't find, and the rule's own log comes before the scheduler's, so a job that died before writing
      its own log (out of memory, never started, etc.) is exactly when the scheduler log doesn't get printed
    - `Logfile PATH ... not found.` (path only)

    Duplicates are dropped, since the same failure is often reported more than once.
    '''
    found = HPCErrors()
    text = [line.rstrip("\n") for line in lines]
    i = 0
    while i < len(text):
        line = text[i].strip()
        i += 1
        if (m := _ERROR_RULE_RE.match(line)):
            if m.group(1) not in found.rules:
                found.rules.append(m.group(1))
        elif (m := _MESSAGE_RE.match(text[i - 1])):
            msg = m.group(1).strip()
            if msg and msg != "None" and msg not in found.messages:
                found.messages.append(msg)
        elif (m := _EXTERNAL_JOBID_RE.match(text[i - 1])):
            if m.group(1) not in found.jobids:
                found.jobids.append(m.group(1))
        elif (m := _LOG_LIST_RE.match(text[i - 1])):
            for path in (x.strip() for x in m.group(1).split(", ")):
                if path and path not in found.listed:
                    found.listed.append(path)
        elif (m := _LOGFILE_HEADER_RE.match(line)):
            # contents sit between two lines of ====, which is the same as "everything until the second one"
            content: list[str] = []
            while i < len(text):
                if is_log_sep(text[i]):
                    i += 1
                    if content:
                        break
                    continue
                content.append(text[i])
                i += 1
            body = "\n".join(content).strip()
            # if the same file shows up more than once, keep the fullest copy
            if len(body) >= len(found.logs.get(m.group(1)) or ""):
                found.logs[m.group(1)] = body
        elif (m := _LOGFILE_EMPTY_RE.match(line)):
            found.logs.setdefault(m.group(1), "")
        elif (m := _LOGFILE_NOTFOUND_PATH_RE.match(line) or _LOGFILE_BINARY_RE.match(line)):
            found.logs.setdefault(m.group(1), None)
    return found


def read_tail(path: str, lines: int = 30, max_bytes: int = 65536) -> str | None:
    '''Last `lines` lines of a text file (reading at most its last `max_bytes` bytes), or None if it can\'t be read'''
    try:
        with open(path, "rb") as f:
            f.seek(0, 2)
            size = f.tell()
            f.seek(max(size - max_bytes, 0))
            data = f.read()
        rows = data.decode("utf-8", errors = "replace").splitlines()
        if size > max_bytes:
            rows = rows[1:]  # the first line is probably cut off
        return "\n".join(rows[-lines:]).strip()
    except OSError:
        return None


class _Pushback:
    '''
    Wraps a single-pass line iterator with one-line lookahead/pushback.
    This is used for processing snakemake output for error handling.
    '''
    __slots__ = ('_it', '_buf')
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
                files = [f.strip() for f in lines[i].split(',')]
                out.extend('  ' + f for f in files)
                i += 1
            continue
        i += 1
    return '\n'.join(out)

def fmt_shell_cmd(text: str, indent: str = "    ") -> str:
    '''
    Reindent a shell block by brace nesting, not by whatever whitespace
    snakemake printed (which drops/adds blank lines and dedents inconsistently).
    '''
    out = []
    depth = 0
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue  # drop snakemake's injected blank lines
        dedent_now = s.startswith('}')
        cur_depth = max(depth - 1, 0) if dedent_now else depth
        out.append(indent * cur_depth + s)
        depth = max(depth + s.count('{') - s.count('}'), 0)
    return '\n'.join(out)


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
        line = next(self.errortext, None)
        if line is None:
            return
        if line.strip().startswith("FileNotFound"):
            if "envs/" in line and ".yaml'" in line:
                self.hp.print('[red]Missing conda environment yaml file:[/][yellow]\n  ' + line.split(':')[-1].replace("'", ''))
            else:
                self.hp.print(line.strip(), style = 'red')
            return

        # ---------- pick out conda-env errors
        if 'Could not create conda' in line:
            for i in self.errortext:
                if 'To search for alternate' in i:
                    break
                if i.strip():
                    if i.lstrip().startswith('-'):
                        self.hp.print(i.rstrip(), soft_wrap = True, width = 2000, style = 'bold red')
                    else:
                        self.hp.print(i.rstrip(), soft_wrap = True, width = 2000, style = 'red')
            return

        if 'singularity image' in line:
            self.hp.print(line.rstrip(), soft_wrap = True, width = 2000, style = 'red')
            for i in self.errortext:
                if i.strip():
                    self.hp.print(i.rstrip(), soft_wrap = True, width = 2000, style = 'red')
            return

        # ---------- pick out snakemake exceptions and missingoutput
        if ('Error' in line or 'Exception' in line or 'Missing input files' in line) and not ('RuleException' in line or 'CalledProcessError' in line):
            self.hp.print(line, highlight=False, soft_wrap = True, end = '', style = 'red')
            for i in self.errortext:
                self.hp.print(i, highlight = False, soft_wrap = True, end = '', style='red')
            return

        if 'but some output files are missing' in line:
            self.missingoutput.append(line)
            for i in self.errortext:
                if 'Shutting down, this might' in i:
                    break
                elif 'but some output files are missing' not in i:
                    self.missingoutput[-1] += i
                else:
                    self.missingoutput.append(i)

        apptainer_store = ''
        for i in self.errortext:
            # ------ it was a command that failed in a container?
            if not i.strip():
                continue

            if "Command ' apptainer  exec --home" in i:
                apptainer_store = i.partition('bash -c')[0] + " bash -c"
                #continue

            hm = _ERROR_RULE_RE.match(i.strip())
            if hm:
                rule = self._parse_rule_block()
                rule.name = hm.group(1)
                if apptainer_store:
                    rule.apptainerprefix = apptainer_store.removeprefix("Command ' ").replace("  ", " ")
                self.rules.append(rule)
                # reset in case of next rule
                apptainer_store = ''

            if i.lstrip().startswith('Logfile '):
                self.errortext.push(i)
                break

            if '(100%) done' in i:
                break

            if 'RuleException' in i:
                #print("ARE WE HITTING THIS EXIT?")
                break
                #return
                #sys.exit(1)

            m = _KEY_RE.match(i)
            if m and m.group(1) == 'message':
                gm = _GROUP_MSG_RE.match(m.group(2).strip())
                if gm:
                    self._skip_group_block()
                    continue

            elif i.startswith('Complete log'):
                break
            elif i.startswith('WorkflowError'):
                break

        if self.missingoutput:
            self.hp.print('\n[bold dim]──── ⚠ Error Reported by Snakemake')
            self.hp.print(fmt_missing_block(self.missingoutput[0]), style = 'red')
        
        elif len(self.rules) > 1:
            grp = ' → '.join(rule.name for rule in self.rules)
            self.hp.rule(
                f'[default bold]Triggering Group[/][yellow bold] {grp}',
                style='yellow'
            )
            for rule in self.rules:
                self.hp.print(f'\n⏺─────── {rule.name}', style = 'bold yellow')
                self.print(rule)

        elif len(self.rules) == 1:
            self.hp.rule(
                f"[default bold]Triggering Rule[/][yellow bold] {self.rules[0].name}",
                style='yellow'
            )
            self.print(self.rules[0])

    def process_hpc(self, directory: str = ".", tail: int = 30, max_logs: int = 4) -> bool:
        '''
        Print what `scrape_hpc_errors` could find in the snakemake output of a workflow run with a scheduler
        plugin, showing the last `tail` lines of up to `max_logs` log files. Log files that were listed but
        that snakemake didn't print are read from disk, with relative paths resolved against `directory` (the
        directory snakemake ran in). Returns False (and prints nothing) if there was nothing useful to show.
        '''
        found = scrape_hpc_errors(self.errortext)
        if not found:
            return False
        # listed logs come first, in order, then any others snakemake printed
        logs: dict[str, str | None] = {path: found.logs.get(path) for path in found.listed}
        for path, body in found.logs.items():
            logs.setdefault(path, body)
        for path, body in logs.items():
            if body is None:
                logs[path] = read_tail(path if os.path.isabs(path) else os.path.join(directory, path), tail)
        found.logs = logs
        self.hp.print('\n[bold dim]──── ⚠ Error Reported by Snakemake [dim](HPC mode, best effort)')
        if found.rules:
            self.hp.print("rule: " + ", ".join(found.rules), style = 'red')
        if found.jobids:
            self.hp.print("scheduler job: " + ", ".join(found.jobids), style = 'red')
        for msg in found.messages:
            self.hp.print(f"message: {msg}", style = 'red', soft_wrap = True, width = 2000, highlight = False)
        for n, (path, body) in enumerate(found.logs.items()):
            if n >= max_logs:
                self.hp.print(f"\n[dim]...and {len(found.logs) - max_logs} more log file(s), see the snakemake log", highlight = False)
                break
            self.hp.print(f"\n──── 🗎 {path}", style = 'bold dim', highlight = False)
            if body is None:
                self.hp.print("(snakemake could not read this file, it may not exist yet or the filesystem hasn't caught up)", style = 'dim')
            elif not body:
                self.hp.print("(empty file)", style = 'red')
            else:
                rows = body.splitlines()
                if len(rows) > tail:
                    self.hp.print(f"[dim]... {len(rows) - tail} earlier line(s) not shown, see the file", highlight = False)
                    rows = rows[-tail:]
                self.hp.print(escape("\n".join(rows)), style = 'red', soft_wrap = True, width = 2000, highlight = False)
        return True

    def _parse_rule_block(self) -> SnakeRule:
        '''Parse one complete Error in rule block.'''
        rule = SnakeRule()
        key = None

        for line in self.errortext:
            if not line.strip():
                continue

            # A logfile belongs to this rule. Stop parsing fields and
            # let _consume_logfile_blocks() handle it.
            if _LOGFILE_HEADER_RE.match(line.strip()):
                self.errortext.push(line)
                break
            # A new rule means this rule has ended.
            if _ERROR_RULE_RE.match(line):
                self.errortext.push(line)
                break
            if line.lstrip().startswith('Complete log'):
                self.errortext.push(line)
                break
            if line.lstrip().startswith('WorkflowError'):
                self.errortext.push(line)
                break

            m = _KEY_RE.match(line)
            if m:
                key, val = m.group(1), m.group(2)

                if key == 'jobid':
                    rule.jobid = int(val)
                elif key == 'input':
                    rule.input = val
                elif key == 'output':
                    rule.output = val
                elif key == 'log':
                    # The logfile block below contains the actual contents,
                    # but retain the path here too.
                    log = val.replace('(check log file(s) for error details)', '').strip()
                    if log:
                        rule.logs.append(log)

                elif key == 'message':
                    rule.message = val
                elif key == 'conda-env':
                    rule.env = val
                elif key == 'resources':
                    rule.resources = [r.strip() for r in val.split(',')]
                elif key == 'shell':
                    rule.cmd = ''

            elif key == 'shell' and line[:1].isspace():
                if line.lstrip().startswith('Logfile'):
                    #self.errortext.push(line)
                    break
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
        '''Consume zero or more Snakemake Logfile blocks.'''
        for line in self.errortext:
            m = _LOGFILE_HEADER_RE.match(line.strip())

            if m:
                logfile = m.group(1)

                content = []

                for j in self.errortext:
                    if is_log_sep(j):
                        if content:
                            break
                        continue
                    content.append(j)

                logtext = escape(
                    re.sub(
                        r'\n{3,}',
                        '\n\n',
                        ''.join(content)
                    ).removeprefix("    ")
                )

                # Purge unnecessary papermill error text
                if "papermill.exceptions.PapermillExecutionError:" in logtext:
                    logtext = logtext.partition('papermill.exceptions.PapermillExecutionError:')[-1]
                    chunks = logtext.split('\n\n')
                    filtered = [
                        c for c in chunks
                        if not c.startswith('File ')
                    ]
                    logtext = '\n\n'.join(filtered)

                # Avoid adding the same logfile twice if the `log:` field
                # was already captured.
                if logfile not in rule.logs:
                    rule.logs.append(logfile)
                rule.log_contents[logfile] = logtext.strip()
                continue

            if _LOGFILE_NOTFOUND_RE.match(line.strip()):
                continue

            self.errortext.push(line)
            break

    def _skip_group_block(self):
        """Skip the incomplete group jobs manifest."""
        for line in self.errortext:
            if _ERROR_RULE_RE.match(line):
                self.errortext.push(line)
                return

    def print(self, rule: SnakeRule):
        'Print a nicely formatted Snakemake rule error'
        self.hp.print("input:", style = 'bold default')
        self.hp.print("  " + rule.input.replace(", ", "\n  "), style = 'red')

        self.hp.print("output:", style = 'bold default')
        self.hp.print("  " + rule.output.replace(", ", "\n  "), style = 'red')

        if rule.logs:
            self.hp.print("log:", style = 'bold default')
            for i in rule.logs:
                _i = i.replace("(check log file(s) for error details)", "").strip()
                self.hp.print("  " + _i, style = 'red')

        if rule.env:
            env_clean = rule.env.split("/")[:-2]
            self.hp.print(
                f"[bold default]conda-env:[/] [red]{'/'.join(env_clean)}[/]"
            )

        if rule.resources:
            for i in rule.resources:
                self.hp.print(i.replace('=', ': '), style = 'red')

        if rule.message != "None":
            self.hp.print("Message:", style = 'bold default')
            self.hp.print(rule.message, style = 'red')

        if rule.cmd:
            
            self.hp.print('')
            self.hp.print("[bold dim]──── ❯ Command Invoked")
            if rule.apptainerprefix:
                self.hp.print(rule.apptainerprefix, end = " '")
            self.hp.shell(rule.cmd, add_quote = bool(rule.apptainerprefix))

        if rule.logs:
            self.hp.print('')

        for i in rule.logs:
            self.hp.print(f"──── 🗎 {i}", style = 'bold dim')
            self.hp.print(rule.log_contents.get(i, "(empty file)"), style = 'red')

