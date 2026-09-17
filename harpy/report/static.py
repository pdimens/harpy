import json
import re
from nbconvert.filters import markdown2html
import shutil
import subprocess
import tempfile
from pathlib import Path

from rich.console import Group

from harpy.common.printing import HarpyPrint

try:
    import yaml
except ImportError:
    yaml = None

FRONTMATTER_RE = re.compile(r"\A---\s*\n(.*?\n)---\s*\n?", re.DOTALL)
DIRECTIVE = re.compile(
    r"^:::\{[^}]*\}[^\n]*\n(?:^:.*\n)*(.*?)\n:::\s*$",
    re.M | re.S,
)
OPEN = re.compile(r"^:::\{(\w+)\}")
REF = re.compile(r"\[\^(\w+)\](?!:)")
DEF_START = re.compile(r"^\[\^(\w+)\]:\s*(.*)$")


def _inline_refs(line: str) -> str:
    return REF.sub(lambda mm: f"<sup>{mm.group(1)}</sup>", line)


def _sanitize(md: str) -> str:
    lines = md.split("\n")
    n_lines = len(lines)
    out, i = [], 0
    while i < n_lines:
        line = lines[i]

        if line.startswith(":::{dropdown}"):
            title = line.lstrip().removeprefix(":::{dropdown}").strip()
            i += 1
            while i < n_lines and re.match(r"^:\w", lines[i]):
                i += 1
            body: list[str] = []
            while i < n_lines and not lines[i].strip().startswith(":::"):
                body.append(lines[i]); i += 1
            i += 1
            inner = markdown2html(_sanitize("\n".join(body))).strip()
            out.append(f"<details>\n<summary>{title}</summary>\n{inner}\n</details>")
            continue

        if line.startswith("```{card}"):
            title = line.lstrip().removeprefix("```{card}").strip()
            i += 1
            body = []
            while i < n_lines and lines[i].strip() != "```":
                body.append(lines[i]); i += 1
            i += 1
            inner = markdown2html(_sanitize("\n".join(body))).strip()
            out.append(inner)
            #out.append(f"<hr>\n{inner}\n")
            continue

        m = OPEN.match(line)
        if m:
            i += 1
            while i < n_lines and re.match(r"^:\w", lines[i]):
                i += 1
            body = []
            while i < n_lines and lines[i].strip() != ":::":
                body.append(lines[i]); i += 1
            i += 1
            inner = markdown2html(_sanitize("\n".join(body))).strip()
            out.append(f"<aside>\n{inner}\n</aside>")
            continue

        m = DEF_START.match(line)
        if m:
            label, rest = m.group(1), m.group(2)
            i += 1
            body = [rest] if rest else []
            while i < n_lines and (lines[i].startswith(("    ", "\t")) or lines[i].strip() == ""):
                if lines[i].strip() == "" and (i + 1 >= n_lines or not lines[i + 1].startswith(("    ", "\t"))):
                    break
                body.append(lines[i].removeprefix("    ").removeprefix("\t"))
                i += 1
            out.append(f"<p><sup>{label}</sup></p>")
            out.append(markdown2html(_sanitize("\n".join(body).strip())).strip())
            continue

        out.append(_inline_refs(line)); i += 1
    return "\n".join(out)

def has_nbconvert():
    hp = HarpyPrint()
    try:
        has = subprocess.run(
            ["jupyter", "nbconvert", "--version"],
            capture_output=True
        ).returncode == 0
    except FileNotFoundError:
        has = False
    if not has:
        _table = hp.table()
        _table.add_column("tool")
        _table.add_column("installation command", style = "green")
        _table.add_row('pip', 'pip install -U nbconvert')
        _table.add_row('conda', 'conda install -c conda-forge nbconvert')
        _table.add_row('pixi', 'pixi add nbconvert')
        hp.error(
            "Missing dependency",
            "jupyter nbconvert is not found on the PATH and is required to proceed.",
            Group("It can be installed using one of these methods:", _table)   
        )

def has_monolith():
    if not shutil.which('monolith'):
        hp = HarpyPrint()
        _table = hp.table()
        _table.add_column("tool")
        _table.add_column("installation", style = "green")
        _table.add_row('cargo', 'cargo install monolith')
        _table.add_row('prebuilt binary', 'add binary to your PATH from https://github.com/Y2Z/monolith/releases')
        hp.error(
            "Missing dependency",
            "Monolith is required to flatten an HTML notebook but was not found on the PATH.",
            Group("Harpy does not provide it, but it can be installed using:", _table)   
        )


class ReportStatic():
    def __init__(self, quiet: bool, static: bool):
        self.quiet: bool = quiet
        self.static: bool = static
        self.hp = HarpyPrint()
        self.hp.console.soft_wrap = True
        self.nbc_log = "ERROR" if quiet else 30
        has_nbconvert()
        if static:
            has_monolith()

    def render_frontmatter_cell(self) -> None:
        """
        If the notebook's first cell is a MyST-style YAML frontmatter block,
        replace it with a plain Markdown header (title/subtitle/date) so
        nbconvert renders something readable instead of the literal '---'
        block. Fields with no display equivalent (e.g. edit_url, which is
        mystmd/site-only) are silently dropped. Any leftover cell content
        after the frontmatter block is preserved below the header.
        """
        cells = self.nb.get("cells", [])
        if not cells:
            return

        first = cells[0]
        source = first.get("source", "")
        if isinstance(source, list):
            source = "".join(source)

        match = FRONTMATTER_RE.match(source)
        if not match:
            return

        if yaml is None:
            self.hp.notice(
                "PyYAML not installed, leaving notebook frontmatter cell as-is "
                "(pip install pyyaml to render it nicely)"
            )
            return

        meta = yaml.safe_load(match.group(1)) or {}

        lines: list[str] = []
        if meta.get("title"):
            lines.append(f"# {meta['title']}")
        if meta.get("subtitle"):
            lines.append(f"### {meta['subtitle']}")
        if meta.get("date"):
            lines.append(f"*{meta['date']}*")

        remainder = source[match.end():].strip()
        if remainder:
            if lines:
                lines.append("")
            lines.append(remainder)

        first["source"] = "\n".join(lines)
        first["cell_type"] = "markdown"
        self.nb['cells'][0] = first


    def sanitize(self, temp_nb_path):
        '''
        Sanitize the mystmd-specific content into plain html that will be properly formatted by nbconvert.
        This includes: dropdowns, footnotes. Serializes to temporary json file `temp_nb_path`.
        '''
        for i, cell in enumerate(self.nb.get("cells", [])):
            if cell["cell_type"] == "markdown":
                src = "".join(cell["source"])
                cell["source"] = _sanitize(src).splitlines(keepends=True)
                self.nb['cells'][i] = cell
        temp_nb_path.write_text(json.dumps(self.nb, indent = 1), encoding="utf-8")

    def run(self, cmd: list[str], **kwargs) -> None:
        if not self.quiet:
            self.hp.log(' '.join(cmd))
        subprocess.run(cmd, check=True, **kwargs)


    def convert(self, notebook: str):
        nb_path: Path = Path(notebook).resolve()
        nb_name = nb_path.stem
        out_path: Path = nb_path.with_name(f"{nb_name}.html")
            
        # Transformed notebook copy lives NEXT TO the original
        # It's cleaned up in `finally` and the original file is never modified.
        tmp_nb_path = nb_path.with_name(f".{nb_name}-tmp.ipynb")
        workdir = Path(tempfile.mkdtemp(prefix="nb2html_"))
        if self.static:
            intermediate_html = workdir / f"{nb_name}.html"
        else:
            intermediate_html = out_path

        try:
            self.nb = json.loads(nb_path.read_text(encoding="utf-8"))
            self.render_frontmatter_cell()
            self.sanitize(tmp_nb_path)
            #tmp_nb_path.write_text(json.dumps(nb_json), encoding="utf-8")
            nbconvert_cmd = [
                "jupyter", "nbconvert", str(tmp_nb_path),
                "--to", "html",
                "--embed-images",
                f"--log-level={self.nbc_log}",
                "--TagRemovePreprocessor.enabled=True",
                "--TagRemovePreprocessor.remove_cell_tags=remove-cell",
                "--TagRemovePreprocessor.remove_input_tags=remove-input",
                "--TagRemovePreprocessor.remove_all_outputs_tags=remove-output",
                "--output", intermediate_html.name,
                "--output-dir", str(intermediate_html.parent),
            ]
            self.run(nbconvert_cmd)

            if not intermediate_html.is_file():
                self.hp.error(
                    "Missing nbconvert output",
                    f"The expected nbconvert output was not found at {intermediate_html}"
                )

            # nbconvert's HTML has no code-split/dynamically-imported JS, so
            # monolith can reliably inline everything it references (CDN
            # vega-embed, MathJax, fonts, etc.) into one working file.
            if self.static:
                monolith_cmd = ["monolith", str(intermediate_html), "-o", str(out_path)]
                if self.quiet:
                    monolith_cmd.append('-q')
                self.run(monolith_cmd)

        finally:
            tmp_nb_path.unlink(missing_ok=True)
            shutil.rmtree(workdir, ignore_errors=True)

        if not self.quiet:
            self.hp.log(f"Done: {out_path}")

