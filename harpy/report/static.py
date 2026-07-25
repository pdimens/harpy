import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from harpy.common.printing import HarpyPrint

try:
    import yaml
except ImportError:
    yaml = None


FRONTMATTER_RE = re.compile(r"\A---\s*\n(.*?\n)---\s*\n?", re.DOTALL)


class ReportStatic():
    def __init__(self, notebook: str, quiet: bool, static: bool):
        self.quiet: bool = quiet
        self.static: bool = static
        self.notebook: str = notebook
        self.hp = HarpyPrint()
        self.hp.console.soft_wrap = True
        self.nbc_log = "ERROR" if quiet else 30

    def render_frontmatter_cell(self, nb: dict) -> None:
        """
        If the notebook's first cell is a MyST-style YAML frontmatter block,
        replace it with a plain Markdown header (title/subtitle/date) so
        nbconvert renders something readable instead of the literal '---'
        block. Fields with no display equivalent (e.g. edit_url, which is
        mystmd/site-only) are silently dropped. Any leftover cell content
        after the frontmatter block is preserved below the header.
        """
        cells = nb.get("cells", [])
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


    def run(self, cmd: list[str], **kwargs) -> None:
        if not self.quiet:
            self.hp.log(' '.join(cmd))
        subprocess.run(cmd, check=True, **kwargs)


    def convert(self):
        nb_path: Path = Path(self.notebook).resolve()
        nb_name = nb_path.stem
        out_path: Path = nb_path.with_name(f"{nb_name}.html")
            

        # Transformed notebook copy lives NEXT TO the original (not in /tmp).
        # It's cleaned up in `finally` and the original file is never modified.
        tmp_nb_path = nb_path.with_name(f".{nb_name}-tmp.ipynb")
        workdir = Path(tempfile.mkdtemp(prefix="nb2html_"))
        if self.static:
            intermediate_html = workdir / f"{nb_name}.html"
        else:
            intermediate_html = out_path

        try:
            nb_json = json.loads(nb_path.read_text(encoding="utf-8"))
            self.render_frontmatter_cell(nb_json)
            tmp_nb_path.write_text(json.dumps(nb_json), encoding="utf-8")
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

