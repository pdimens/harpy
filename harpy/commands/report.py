import copy
import glob
import os
from pathlib import Path
from time import sleep
import re
import rich_click as click
from rich.live import Live
from rich.panel import Panel
import shutil
import subprocess
import sys
import uuid
import yaml
from harpy.common.file_ops import fetch_template
from harpy.common.printing import print_error, print_notice, CONSOLE

class ReportRender():
    def __init__(self, root: str = ""):
        self.root = root
        self.configfile = os.path.join(root, "myst.yml")
        self.filetree: dict = {}
        self.filechanges: bool = False

        if not os.path.exists(self.configfile):
            _yml = myst_yaml()
            self.scan_for_reports()
            with open(self.configfile, "w") as yml:
                yaml.dump(_yml, yml, default_flow_style= False, sort_keys=False, width=float('inf'))
        else:
            with open(self.configfile, "r") as yml:
                _yml = yaml.full_load(yml)

        if not isinstance(_yml, dict) or "project" not in _yml:
            print_error(
                "invalid MyST configuration",
                "The [blue]myst.yml[/] file that was found is improperly formatted.",
                "If you manually edited this file, please try to remake it using [blue]harpy report init[/], otherwise check that the file is in proper YAML format and has a [green]project[/] section."
            )
        if "toc" in _yml["project"]:
            self.toc_to_filetree(_yml["project"].pop("toc"))
        self.config: dict = _yml

    def update_yaml(self):
        """
        Convert self.filelist to a MyST-formatted table of contents and overwrite `self.configfile` with the updated table of contents
        """
        if not self.filechanges:
            return
        self.config["project"]["toc"] = self.filetree_to_toc()
        with open(self.configfile, "w") as yml:
            yaml.dump(self.config, yml, default_flow_style= False, sort_keys=False, width=float('inf'))
        del self.config["project"]["toc"]

    def scan_for_reports(self):
        """
        Recursively search `self.root` for files ending in `.ipynb`, filtering out those found in the _build/ directory
        and converts that list of file paths into a nested dictionary tree structure stored as `self.filetree`.
        """
        original = copy.deepcopy(self.filetree)
        _ipynb = set(i for i in glob.iglob("**/*.ipynb", root_dir = self.root, recursive = True) if "_build" not in i and "workflow/report" not in i)
        for path in _ipynb:
            # ignore the 'reports' folder name when building the tree
            parts = [i for i in Path(path).parts if i != "reports"]

            if len(parts) == 1:
                # Root level file
                if '__root__' not in self.filetree:
                    self.filetree['__root__'] = []
                self.filetree['__root__'].append(parts[0])
            else:
                # Navigate to parent directory, creating dicts as needed
                current = self.filetree
                for part in parts[:-2]:
                    if part not in current:
                        current[part] = {}
                    current = current[part]

                # Add the file to its parent directory
                parent_dir = parts[-2]

                if parent_dir not in current:
                    current[parent_dir] = set()

                # If it's already a set, add to it
                if isinstance(current[parent_dir], set):
                    current[parent_dir].add(path)
                # If it's a dict, it has subdirs, so we need to keep it as dict
                # (this handles cases where a dir has both files and subdirs)
        if original != self.filetree:
            self.clean_filetree()
            self.filechanges = True
        else:
            self.filechanges = False

    def clean_filetree(self):
        """
        Scans the filetree for nonexistant files and removes them, then recusively cleans up keys that terminate
        in an empty value
        """
        def remove_nonexistent_files(d):
            if not isinstance(d, dict):
                return d
            cleaned = {}
            for key, value in d.items():
                if isinstance(value, dict):
                    # Recursively process nested dictionaries
                    cleaned_value = remove_nonexistent_files(value)
                    # Only keep if the cleaned dict is not empty
                    if cleaned_value:
                        cleaned[key] = cleaned_value
                elif isinstance(value, set):
                    # Filter out non-existent files from the list
                    existing_files = set(f for f in value if os.path.exists(f))
                    # Only keep the key if there are existing files
                    if existing_files:
                        cleaned[key] = existing_files
                else:
                    # For non-dict, non-list values, keep as is (unless explicitly empty)
                    if value: # not in (None, '', [], {}, (), set()):
                        cleaned[key] = value
            return cleaned
        self.filetree = remove_nonexistent_files(self.filetree)

    def filetree_to_toc(self):
        """
        Recursively parses `self.filetree` to format it as a MyST table of contents and return it
        """
        def recursive_transform(tree):
            """
            Transform a nested dict into a format where dict values become 
            {'title': key, 'children': transformed_value}.
            Terminal lists of filenames become lists of {'file': filename} dicts.
            
            Returns a transformed structure with title/children format
            """
            if isinstance(tree, dict):
                result = []
                if tree.get("root", []):
                    result += [{'file': filename} for filename in tree.get("root", [])]

                for key, value in tree.items():
                    if key == 'root':
                        continue    
                    if isinstance(value, dict):
                        # Recursively transform the nested dict
                        result.append({
                            'title': key,
                            'children': recursive_transform(value)
                        })
                    else:
                        # Value is a list (leaf node) - convert filenames to {'file': filename}
                        result.append({
                            'title': key,
                            'children': [{'file': filename} for filename in value]
                        })
                return result
            else:
                # If it's not a dict (shouldn't happen at top level), return as-is
                return tree

        return recursive_transform(self.filetree)

    def toc_to_filetree(self, toc):
        """
        Recursively parses a `toc` (as interpreted in myst.yml) to format it as a nested `dict` stored in `self.filetree`
        """
        def recursive_transform(transformed):
            if isinstance(transformed, list):
                result = {}
                root_files = []

                for item in transformed:
                    # Handle case where item is just {"file": filename} at top level
                    if 'file' in item and 'title' not in item:
                        root_files.append(item['file'])
                        continue

                    title = item['title']
                    children = item['children']

                    # If children is a list of dicts with 'title' key, recurse
                    if isinstance(children, list) and children and isinstance(children[0], dict):
                        if 'title' in children[0]:
                            # Nested structure
                            result[title] = recursive_transform(children)
                        elif 'file' in children[0]:
                            # Terminal list of files - extract filenames
                            result[title] = [child['file'] for child in children]
                    else:
                        # Children is some other structure (shouldn't happen with our format)
                        result[title] = children

                # Add root files if any were found
                if root_files:
                    result['root'] = root_files

                return result
            else:
                # If it's not a list, return as-is
                return transformed

        # merge the dicts
        self.filetree = recursive_transform(toc)

def rand_id() -> str:
    """
    Returns a random 32-digit hyphenated uuid with lengths 8-4-4-4-12,
    e.g. `c7771464-3fd4-42c9-b6d1-45a1c5c656e1`
    """
    uid = uuid.uuid4().hex
    indices = [0, 8, 12, 16, 20, 32]
    return "-".join(uid[i:j] for i,j in zip(indices, indices[1:]))

def myst_yaml() -> dict:
    """
    Returns a dict with the contents of a harpy-configured `mysy.yml` file
    """
    git_url = subprocess.run("git remote get-url origin".split(), text = True, stdout = subprocess.PIPE).stdout.strip().removesuffix(".git")
    d = {
        "version" : 1,
        "site" : {
            "template": "book-theme",
            "options" : {
                "favicon" : "PLACEHOLDER",
                "logo" : "PLACEHOLDER"
            },
        },
        "project" : {
            "id": rand_id(),
            "github" : git_url if git_url else {},
            "title" : "Harpy Reports",
            "description" : "The reports produced by Harpy, aggregated into a navigable website using MyST.",
            "toc": [{"file" : ".report/index.md"}]
        }
    }
    return d

@click.group(context_settings={"help_option_names" : []})
def report():
    """
    Setup and render harpy reports

    Harpy reports are provided as Jupyter Notebooks. The subsequent commands
    will configure your project directory to make these much nicer, enabling you
    to render them as one navigable and interactive local website.
    """

@click.command(context_settings={"help_option_names" : ['-h', '--help']}, epilog = "Documentation: https://pdimens.github.io/harpy/reports/")
def init():
    """
    Configure the project for advanced report features
    
    Using MyST, Harpy reports can be built into and rendered as
    an interactive website, which first requires a bit of configuration.
    """
    notices = []
    if not shutil.which("git"):
        notices.append("- The [green]git[/] software was not found to identify the root directory of this project, so it was assumed this is the root directory.")
        git_dir = ""
    else:
        git_dir = subprocess.run("git rev-parse --show-toplevel".split(), text = True, stdout = subprocess.PIPE).stdout.strip()
    #TODO pull any specific configurations here too
    fetch_template("report_index.md", os.path.join(git_dir, ".report", "index.md"))
    with open(os.path.join(git_dir, "myst.yml"), "w") as yml:
        yaml.dump(myst_yaml(), yml, default_flow_style= False, sort_keys=False, width=float('inf'))

    if not shutil.which("myst"):
        notices.append("- The [green]MyST[/] software is required for a locally-served report website, however it was not found in this environment.")
    if notices:
        print_notice("\n".join(notices))

@click.command(context_settings={"help_option_names" : ['--help']}, epilog = "Documentation: https://pdimens.github.io/harpy/reports/")
@click.option('-h', '--headless', is_flag = True, help = 'Run the server in headless mode, with only the content server started')
@click.option('-p', '--port', type = int, help = 'Run the application server from the specified port number')
@click.option('-r', '--refresh', type = click.IntRange(min = 0, max_open=True), show_default = True, default = 30, help = 'Refresh interval, in seconds')
@click.option('-s', '--server-port', type = int, help = 'Run the content server from the specified port number')
@click.argument('directory', required=False, type = click.Path(exists = True, file_okay = False, readable = True), nargs = 1)
def render(directory, headless, port, server_port, refresh):
    """
    Render ipynb reports as a local website

    Using MyST, all the `.ipynb` reports within Harpy-generated
    directories will be aggregated and rendered into a locally-served
    website for you to review them from a single access point. This command
    is expected to be executed within a git version-controlled directory, where
    Harpy can identify the root directory of the project, otherwise provide the
    path to a directory for Harpy to recursively scan the `.ipynb` reports. 
    """
    cmd = ["myst", "start"]
    if headless:
        cmd.append("--headless")
    if port:
        cmd += ["--port", f"{port}"]
    if server_port:
        cmd += ["--server-port", f"{server_port}"]

    tracker = ReportRender(directory if directory else "")
    tracker.scan_for_reports()
    tracker.filetree_to_toc()
    tracker.update_yaml()
    URL = ""
    myst_error = ""
    try:
        panel = Panel("Starting the MyST live-server[dim]...", border_style = "medium_purple4", title = "[default bold]Harpy report", subtitle= "[default]Terminate it with[/] [bold yellow]ctrl+c[/]")
        with subprocess.Popen(cmd, cwd = directory, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True) as serve, Live(panel, console = CONSOLE, auto_refresh = False, transient = True) as live:
            while not URL:
                if serve.poll():
                    myst_error += "\n".join(serve.stderr.readlines())
                    raise ValueError
                _myst_output = serve.stdout.readline()
                _url = re.findall(r'http://\S+\s', _myst_output)
                if _url:
                    URL += _url[0].strip()
            panel.renderable = f"MyST live-server started: [blue bold]{URL}[/]"
            live.refresh()
            while True:
                if serve.poll():
                    raise ValueError
                if refresh > 0:
                    tracker.scan_for_reports()
                    tracker.update_yaml()
                    sleep(refresh)
                else:
                    sleep(999999)
    except KeyboardInterrupt:
        # clear the top part of the panel
        for _ in range(1):
            CONSOLE.file.write("\033[F\033[K")
        CONSOLE.file.flush()
    except ValueError:
        print_error(
            "MyST server error",
            f"The [blue]myst start[/] command exited and reported this error:\n[yellow]{myst_error.strip()}[/]"
        )

report.add_command(init)
report.add_command(render)
