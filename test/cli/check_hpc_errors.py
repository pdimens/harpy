"""
Checks for the best-effort error reporting of workflows run with a scheduler (HPC) plugin. Run with plain python:

    python test/cli/check_hpc_errors.py

The snakemake output below follows the formats in snakemake's own error reporting and in the slurm, lsf,
googlebatch, and cluster-generic executor plugins. Real output from `cluster-generic` was used as a reference.
"""

import io
import os
import sys
import tempfile

from harpy.common.errorparsing import ErrorHandler, read_tail, scrape_hpc_errors

FAILED = []

def check(condition: bool, message: str) -> None:
    print(f"{'ok  ' if condition else 'FAIL'} {message}")
    if not condition:
        FAILED.append(message)

def printed(lines, directory=".", **kwargs) -> tuple[bool, str]:
    """run the handler and capture what it prints"""
    handler = ErrorHandler(lines)
    handler.hp.console.file = io.StringIO()
    handler.hp.console.width = 120
    result = handler.process_hpc(directory, **kwargs)
    return result, handler.hp.console.file.getvalue()

SEP = "=" * 60

def block(path: str, *content: str) -> list[str]:
    """how snakemake prints a log file's contents"""
    return [f"Logfile {path}:\n", SEP + "\n", *[c + "\n" for c in content], SEP + "\n"]

# ---- slurm: both logs exist and snakemake printed both
SLURM_LOG = ".snakemake/slurm_logs/rule_bwa_map/s1/12345.log"
slurm = [
    "[Wed Sep 30 19:45:44 2026]\n",
    "Error in rule bwa_map:\n",
    "    message: SLURM-job '12345' failed, SLURM status is: 'FAILED'. For further error details see the cluster/cloud log and the log files of the involved rule(s).\n",
    "    jobid: 4\n",
    "    input: genome.fa, s1.fq\n",
    "    output: bam/s1.bam\n",
    f"    log: logs/bwa/s1.log, {SLURM_LOG} (check log file(s) for error details)\n",
    "    shell:\n        bwa mem genome.fa s1.fq > bam/s1.bam\n        (command exited with non-zero exit code)\n",
    "    external_jobid: 12345\n",
    *block("logs/bwa/s1.log", "[M::bwa_idx_load] reading index", "[E::main] fail to open file"),
    *block(SLURM_LOG, "slurmstepd: error: Detected 1 oom_kill event in StepId=12345.batch."),
    "Shutting down, this might take some time.\n",
]
found = scrape_hpc_errors(slurm)
check(found.rules == ["bwa_map"], "slurm: rule name")
check(found.jobids == ["12345"], "slurm: scheduler job id")
check(len(found.messages) == 1 and "SLURM status is: 'FAILED'" in found.messages[0], "slurm: message")
check(list(found.logs) == ["logs/bwa/s1.log", SLURM_LOG], "slurm: both log files")
check("oom_kill" in found.logs[SLURM_LOG], "slurm: scheduler log contents")
ok, out = printed(slurm)
check(ok and "oom_kill" in out and "fail to open file" in out and "12345" in out, "slurm: printed output has the useful text")

# ---- the first line (the "Error in rule" header) is consumed by harpy's monitor and not always available
no_header = slurm[2:]
found = scrape_hpc_errors(no_header)
check(found.rules == [] and "oom_kill" in found.logs[SLURM_LOG] and found.messages, "header line missing: everything else is still found")
found = scrape_hpc_errors(["Error in rule bwa_map:\n"] + no_header)
check(found.rules == ["bwa_map"], "header line prepended (as workflow.py does): rule name is found")

# ---- the rule log was never written, so snakemake stops at "not found" and never prints the scheduler log
with tempfile.TemporaryDirectory() as outdir:
    os.makedirs(os.path.join(outdir, ".snakemake/slurm_logs/rule_bwa_map/s1"))
    with open(os.path.join(outdir, SLURM_LOG), "w") as f:
        f.write("slurmstepd: error: *** JOB 12345 ON node1 CANCELLED DUE TO TIME LIMIT ***\n")
    missing_rule_log = [l for l in slurm if not l.startswith("Logfile") and set(l.strip()) != {"="}
                        and "oom_kill" not in l and "[M::" not in l and "[E::" not in l]
    missing_rule_log.insert(-1, "Logfile logs/bwa/s1.log not found.\n")
    check(SLURM_LOG not in scrape_hpc_errors(missing_rule_log).logs or scrape_hpc_errors(missing_rule_log).logs[SLURM_LOG] is None,
          "rule log missing: snakemake itself didn't print the scheduler log")
    ok, out = printed(missing_rule_log, outdir)
    check(ok and "CANCELLED DUE TO TIME LIMIT" in out, "rule log missing: the scheduler log is read from disk")
    check("could not read" in out, "rule log missing: says the missing rule log couldn't be read")

    # relative paths are resolved against the directory snakemake ran in, not the cwd
    os.makedirs(os.path.join(outdir, "logs"))
    with open(os.path.join(outdir, "logs/x.log"), "w") as f:
        f.write("a relative log\n")
    ok, out = printed(["    log: logs/x.log (check log file(s) for error details)\n"], outdir)
    check(ok and "a relative log" in out, "relative log paths are resolved against the snakemake directory")

# ---- lsf: same shape
LSF_LOG = "/home/u/proj/.snakemake/lsf_logs/rule_qc/s2/9f8e.log"
lsf = [
    "Error in rule qc:\n",
    "    message: LSF job 778 failed, job status: EXIT\n",
    "    jobid: 2\n",
    f"    log: logs/qc/s2.log, {LSF_LOG} (check log file(s) for error details)\n",
    *block(LSF_LOG, "TERM_MEMLIMIT: job killed after reaching LSF memory usage limit."),
]
found = scrape_hpc_errors(lsf)
check(found.rules == ["qc"] and LSF_LOG in found.logs and "TERM_MEMLIMIT" in found.logs[LSF_LOG], "lsf: rule, scheduler log, and contents")

# ---- google batch: only the plugin's log (fetched from cloud logging) plus the rule log
GB_LOG = ".snakemake/googlebatch_logs/rule_align/s3.log"
gb = [
    "Error in rule align:\n",
    "    message: Google Batch job 'projects/p/locations/us-central1/jobs/j-1' failed. \n",
    "    jobid: 1\n",
    f"    log: {GB_LOG} (check log file(s) for error details)\n",
    *block(GB_LOG, "Task failed: exit code 137"),
]
found = scrape_hpc_errors(gb)
check("Task failed" in found.logs[GB_LOG] and "failed." in found.messages[0], "googlebatch: plugin log and message")

# ---- cluster-generic, abbreviated from real snakemake 9.27 output: nested job output is interleaved with the
# ---- head node's blocks, so the same failure is reported twice with different job ids
generic = [
    "RuleException:\n",
    "CalledProcessError in file \"/p/Snakefile\", line 8:\n",
    "Command 'set -euo pipefail; ls /nonexistent' returned non-zero exit status 2.\n",
    "Error in rule work:\n",
    "    message: None\n", "    jobid: 0\n", "    output: out/a.txt\n",
    "    log: logs/work.a.log (check log file(s) for error details)\n",
    "Exiting because a job execution failed. Look above for error messages\n",
    "Error in rule work:\n",
    "    message: Error submitting jobscript (exit code 1):\n",
    "For further error details see the cluster/cloud log and the log files of the involved rule(s).\n",
    "    jobid: 1\n", "    output: out/a.txt\n",
    "    log: logs/work.a.log (check log file(s) for error details)\n",
    *block("logs/work.a.log", "starting sample a", "ls: cannot access '/nonexistent': No such file or directory"),
]
found = scrape_hpc_errors(generic)
check(found.rules == ["work"], "cluster-generic: the repeated failure is reported once")
check(found.messages == ["Error submitting jobscript (exit code 1):"], "cluster-generic: 'message: None' is ignored and the message isn't duplicated")
check(list(found.logs) == ["logs/work.a.log"] and "nonexistent" in found.logs["logs/work.a.log"], "cluster-generic: rule log contents")

# ---- snakemake's other log outcomes
found = scrape_hpc_errors(["Logfile logs/e.log: empty file\n", "Logfile logs/b.bin is not a text file.\n"])
check(found.logs == {"logs/e.log": "", "logs/b.bin": None}, "empty and non-text log files are recognized")

# ---- long logs are cut to the last lines
long_log = block(SLURM_LOG, *[f"line {n}" for n in range(1, 101)])
ok, out = printed(long_log, tail=10)
check("line 100" in out and "line 91" in out and "line 90\n" not in out and "90 earlier line(s) not shown" in out, "long logs show only the last lines, and say so")
with tempfile.TemporaryDirectory() as d:
    path = os.path.join(d, "big.log")
    with open(path, "w") as f:
        f.write("\n".join(f"row {n}" for n in range(200000)) + "\n")
    tail = read_tail(path, 3)
    check(tail == "row 199997\nrow 199998\nrow 199999", "read_tail reads the end of a large file without loading all of it")
check(read_tail("/nonexistent/file.log") is None, "read_tail returns None for unreadable files")

# ---- markup in logs must not be interpreted
ok, out = printed(block("logs/m.log", "[bold red]not markup[/] and [/ broken"))
check(ok and "[bold red]not markup[/] and [/ broken" in out, "rich markup in log contents is printed literally")

# ---- nothing useful: the caller prints its own message
ok, out = printed(["Building DAG of jobs...\n", "Shutting down, this might take some time.\n", "Exiting because a job execution failed.\n"])
check(not ok and out == "", "nothing found: returns False and prints nothing")

if FAILED:
    print(f"\n{len(FAILED)} check(s) failed")
    sys.exit(1)
print("\nall checks passed")
