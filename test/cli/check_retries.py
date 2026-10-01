"""
Checks that LaunchSnakemake tolerates snakemake retrying failed jobs (--retries, `retries:` in HPC profiles) and only
reports a failure once snakemake has really given up. Run with plain python:

    python test/cli/check_retries.py

Snakemake is replaced by a process that replays recorded output, so snakemake isn't needed. The files in
test/cli/fixtures were recorded from snakemake 9.27 with a job that fails once and then succeeds, a job that fails every
attempt (--retries 2), and a job that fails without retries. The other cases are written out here, following the
formats of snakemake's scheduler and of the scheduler plugins.
"""

import io
import os
import sys
import tempfile
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

from harpy.common.errorparsing import ErrorHandler
from harpy.common.launch import LaunchSnakemake
from harpy.common.printing import HarpyPrint

FIXTURES = Path(__file__).resolve().parent / "fixtures"
FAILED = []

def check(condition: bool, message: str) -> None:
    print(f"{'ok  ' if condition else 'FAIL'} {message}")
    if not condition:
        FAILED.append(message)

REPLAY = '''
import sys
sys.stderr.write(open(sys.argv[1]).read())
sys.stderr.flush()
sys.exit(int(sys.argv[2]))
'''

def run(tmp: str, lines: list[str], exit_code: int, quiet: int = 2) -> LaunchSnakemake:
    """Run LaunchSnakemake on a process that prints `lines` to stderr and exits with `exit_code`"""
    replay, output = os.path.join(tmp, "replay.py"), os.path.join(tmp, "stderr.txt")
    Path(replay).write_text(REPLAY)
    Path(output).write_text("".join(lines))
    with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()):
        return LaunchSnakemake(f"{sys.executable} {replay} {output} {exit_code}", tmp, quiet, HarpyPrint())

def fixture(name: str) -> list[str]:
    return (FIXTURES / name).read_text().splitlines(keepends=True)

def blocks(sm: LaunchSnakemake) -> int:
    """number of "Error in rule" blocks passed on (the line that triggered the error is `output`, the rest `errorlog`)"""
    return sum(1 for line in [sm.output] + sm.errorlog if line.strip().startswith("Error in rule"))

# the start of every snakemake run, up to the first jobs being selected, from the recorded output
PRELUDE = []
for line in fixture("retry_then_success.err"):
    PRELUDE.append(line)
    if line.startswith("Select jobs to execute"):
        break

with tempfile.TemporaryDirectory() as tmp:
    # ---- recorded: a job fails once and succeeds on its retry
    lines = fixture("retry_then_success.err")
    sm = run(tmp, lines, 0)
    check(sm.exitcode == 0, "recorded: a job that fails and then succeeds on retry is not a failure")
    check(sm.retries == {1: 1}, "recorded: the retry is noted")
    check(not sm.errorlog, "recorded: nothing is left over as an error")
    sm = run(tmp, lines, 0, quiet=1)
    done = {name: (sm.progress.tasks[task].completed, sm.progress.tasks[task].total) for name, task in sm.task_ids.items()}
    check(done["flaky"] == (1, 1) and done["steady"] == (1, 1), "recorded: the progress bar goes on through the retry and completes")

    # ---- recorded: every attempt fails. The error of the last attempt is what is reported, once.
    sm = run(tmp, fixture("retries_exhausted.err"), 1)
    check(sm.exitcode == 3, "recorded: a job that fails every attempt is a failure")
    check(sm.retries == {1: 2}, "recorded: both restarts are noted")
    check(blocks(sm) == 1, "recorded: only the error of the last attempt is passed on")
    check(sm.output.strip() in ("RuleException:", "Error in rule never:"), "recorded: the first line of the error is the output, as before")
    check(any("Logfile logs/never.log" in line for line in sm.errorlog), "recorded: the log of the failure is included")
    handler = ErrorHandler(sm.errorlog)
    handler.hp.console.file = io.StringIO()
    handler.process()
    text = handler.hp.console.file.getvalue()
    check("Triggering Rule never" in text and "Triggering Group" not in text, "recorded: the parser reports one failing rule, not a chain of attempts")

    # ---- recorded: no retries at all, which is how it was before
    sm = run(tmp, fixture("no_retries_failure.err"), 1)
    check(sm.exitcode == 3 and blocks(sm) == 1 and not sm.retries, "recorded: a failure without retries is a failure straight away")

    # ---- scheduler plugins: a failed job, then the restart (the first error line is the "Error in rule" header itself)
    def error_block(jobid: int, rule: str = "flaky") -> list[str]:
        return [
            "[Thu Oct  1 12:00:00 2026]\n",
            f"Error in rule {rule}:\n",
            "    message: SLURM-job '12345' failed, SLURM status is: 'FAILED'. For further error details see the cluster/cloud log.\n",
            f"    jobid: {jobid}\n",
            f"    output: out/{rule}.txt\n",
            f"    log: logs/{rule}.log, .snakemake/slurm_logs/rule_{rule}/12345.log (check log file(s) for error details)\n",
            "Logfile .snakemake/slurm_logs/rule_flaky/12345.log:\n", "=" * 40 + "\n", "slurmstepd: error: Detected 1 oom_kill event\n", "=" * 40 + "\n",
        ]
    started = lambda jobid, rule="flaky": ["[Thu Oct  1 12:00:00 2026]\n", f"rule {rule}:\n", f"    output: out/{rule}.txt\n", f"    jobid: {jobid}\n"]
    finished = ["Finished jobid: 1 (Rule: flaky)\n", "1 of 3 steps (33%) done\n", "Finished jobid: 2 (Rule: steady)\n", "2 of 3 steps (67%) done\n",
                "Finished jobid: 0 (Rule: all)\n", "3 of 3 steps (100%) done\n", "Complete log(s): /work/.snakemake/log/x.snakemake.log\n"]
    sm = run(tmp, PRELUDE + started(1) + error_block(1) + ["Trying to restart job 1.\n", "Select jobs to execute...\n", "Execute 1 jobs...\n"] + started(1) + finished, 0)
    check(sm.exitcode == 0 and sm.retries == {1: 1}, "scheduler: a failed job that is restarted and succeeds is not a failure")

    sm = run(tmp, PRELUDE + started(1) + error_block(1) + ["Trying to restart job 1.\n"] + started(1) + error_block(1)
             + ["Finished jobid: 2 (Rule: steady)\n", "Shutting down, this might take some time.\n", "Exiting because a job execution failed. Look above for error messages\n"], 1)
    check(sm.exitcode == 3 and blocks(sm) == 1 and sm.retries == {1: 1}, "scheduler: out of retries is a failure, reported once")
    check(sm.output.strip() == "Error in rule flaky:", "scheduler: the output is the 'Error in rule' header, which the HPC parser gets")

    # ---- one job fails for good while another is restarted: the restart mustn't hide the failure
    sm = run(tmp, PRELUDE + started(2, "steady") + error_block(2, "steady") + error_block(1) + ["Trying to restart job 1.\n"]
             + ["Shutting down, this might take some time.\n", "Exiting because a job execution failed. Look above for error messages\n"], 1)
    check(sm.exitcode == 3, "concurrent: a restart of one job doesn't hide another job's failure")
    check(any(line.strip() == "jobid: 2" for line in sm.errorlog), "concurrent: the error of the job that failed for good is kept")

    # ---- a restart that can't be matched to a job (e.g. a group) must not turn into a failure by itself
    group = ["Error in group g1:\n", "    message: None\n", "    jobs:\n", "        rule a:\n", "            jobid: 5\n"]
    sm = run(tmp, PRELUDE + started(5) + group + ["Trying to restart job 7.\n"] + started(5) + finished[:-2], 0)
    check(sm.exitcode == 0, "group: a restart that isn't matched to a job isn't a failure when snakemake goes on and exits fine")

    # ---- text that only looks like an error, in a run that goes fine
    sm = run(tmp, PRELUDE + started(1) + ["Warning: ErrorDocument is missing from the config\n"] + finished[:-2], 0)
    check(sm.exitcode == 0 and not sm.errorlog, "text with 'Error' in a run that exits fine isn't a failure")

    # ---- no verdict from snakemake, it just stops
    sm = run(tmp, PRELUDE + started(1) + ["Traceback (most recent call last):\n", "KeyboardInterrupt\n", "Error: it died\n"], 1)
    check(sm.exitcode == 3 and any("it died" in line for line in [sm.output] + sm.errorlog), "an error followed by snakemake exiting without the usual markers is a failure")

if FAILED:
    print(f"\n{len(FAILED)} check(s) failed")
    sys.exit(1)
print("\nall checks passed")
