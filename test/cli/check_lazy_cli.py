"""
Checks for the lazily-loaded CLI and shell completion. Run with plain python (no pytest needed):

    python test/cli/check_lazy_cli.py

- the lazy manifests (`COMMANDS`) match the real commands
- importing the entry points doesn't import the command modules or heavy dependencies
- shell completion scripts can be generated and complete what they should
"""

import inspect
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import click
from click.shell_completion import ShellComplete

ROOT = Path(__file__).resolve().parents[2]
FAILED = []

def check(condition: bool, message: str) -> None:
    print(f"{'ok  ' if condition else 'FAIL'} {message}")
    if not condition:
        FAILED.append(message)

def first_paragraph(command: click.Command) -> str:
    return inspect.cleandoc(command.help or "").split("\n\n")[0].replace("\n", " ").strip()

# --- laziness: run in a fresh interpreter so nothing else has been imported
def imported_after(statement: str) -> set[str]:
    code = f"import sys; {statement}; print('\\n'.join(sys.modules))"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    return set(out.stdout.split())

BANNED = ("pysam", "yaml", "nbconvert", "polars", "numpy", "harpy.commands", "harpy.validation", "harpy.common.workflow")
for entrypoint in ("harpy.__main__", "harpy.utils"):
    mods = imported_after(f"import {entrypoint}")
    leaked = sorted(m for m in mods if m.split(".")[0] in BANNED or any(m.startswith(b + ".") or m == b for b in BANNED))
    check(not leaked, f"importing {entrypoint} imports no command modules or heavy dependencies" + (f" (leaked: {leaked[:5]})" if leaked else ""))

# --- manifests match the real commands
from harpy.__main__ import cli, COMMANDS
import harpy.utils

for group, manifest, has_help in ((cli, COMMANDS, True), (harpy.utils.cli, harpy.utils.COMMANDS, False)):
    ctx = click.Context(group)
    for name, spec in manifest.items():
        real = group.resolve_command(ctx, [name])[1]
        check(real.name == name, f"[{group.name}] '{name}' resolves to a command named '{real.name}'")
        if has_help:
            check(spec.help == first_paragraph(real), f"[{group.name}] '{name}' help text matches its docstring")
            check(spec.hidden == real.hidden, f"[{group.name}] '{name}' hidden flag matches")
    # everything defined in the module is in the manifest and nothing extra
    check(sorted(group.list_commands(ctx)) == sorted(manifest), f"[{group.name}] list_commands matches the manifest")

# --- shell completion
def complete(args: list[str], incomplete: str) -> list[tuple[str, str]]:
    comp = ShellComplete(cli, {}, "harpy", "_HARPY_COMPLETE")
    return [(i.type, i.value) for i in comp.get_completions(args, incomplete)]

values = [v for _, v in complete([], "")]
check("align" in values and "qc" in values, "completion offers subcommands")
check("completion" not in values and "containerize" not in values, "completion doesn't offer hidden subcommands")
check([v for _, v in complete([], "al")] == ["align"], "completion narrows subcommands by prefix")
check({"bwa", "strobe", "minimap"} <= {v for _, v in complete(["align"], "")}, "completion descends into nested groups")
check("--threads" in {v for _, v in complete(["align", "bwa"], "--th")}, "completion offers options")
check(("file", "") in complete(["align", "bwa"], ""), "completion of file-like arguments defers to the shell")

for shell in ("bash", "zsh", "fish"):
    out = subprocess.run([sys.executable, "-m", "harpy", "completion", shell], capture_output=True, text=True)
    check(out.returncode == 0 and "_HARPY_COMPLETE" in out.stdout, f"`harpy completion {shell}` prints a completion script")

# --- the conda recipe only has build.sh, so it carries its own copy of the activation hooks
build_sh = (ROOT / "resources" / "build.sh").read_text()
for marker, resource in (("EOF_COMPLETION_SH", "shell_completion.sh"), ("EOF_COMPLETION_FISH", "shell_completion.fish")):
    match = re.search(rf"<<'{marker}'\n(.*?)\n{marker}\n", build_sh, re.S)
    check(bool(match) and match.group(1).strip() == (ROOT / "resources" / resource).read_text().strip(),
          f"resources/build.sh embeds an identical copy of resources/{resource}")

# --- the activation hook must never export FPATH: an FPATH in the environment replaces zsh's whole
# function path (and pixi re-applies it after ~/.zshrc), which breaks oh-my-zsh, prompts, etc.
hook = ROOT / "resources" / "shell_completion.sh"
env = {k: v for k, v in os.environ.items() if k not in ("FPATH", "XDG_DATA_DIRS")}
with tempfile.TemporaryDirectory() as prefix:
    (Path(prefix) / "share").mkdir()
    probe = subprocess.run(["bash", "-c", f'CONDA_PREFIX={prefix} . {hook}; echo "FPATH=${{FPATH-unset}}"; echo "XDG=${{XDG_DATA_DIRS-unset}}"'],
                           capture_output=True, text=True, env=env)
check("FPATH=unset" in probe.stdout, "activation hook does not export FPATH")
check(f"XDG=/usr/local/share:/usr/share:{prefix}/share" in probe.stdout, "activation hook adds the environment to XDG_DATA_DIRS")

if FAILED:
    print(f"\n{len(FAILED)} check(s) failed")
    sys.exit(1)
print("\nall checks passed")
