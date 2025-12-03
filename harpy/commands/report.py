import glob
import os
from pathlib import Path
from time import sleep
import re
import rich_click as click
from rich import print as rprint
import shutil
import subprocess
import sys
import uuid
import yaml
from harpy.common.file_ops import fetch_template
from harpy.common.printing import print_error, print_notice

class ReportRender():
    def __init__(self, root: str = ""):
        self.root = root
        self.configfile = os.path.join(root, "myst.yml")
        self.filetree: dict = {}
        self.filechanges: bool = False
        #self.toc = []
        if not os.path.exists(self.configfile):
            _yml = myst_yaml()
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

    def overwrite_yaml(self):
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
        tree = self.filetree
        _ipynb = set(i for i in glob.iglob("**/*.ipynb", root_dir = self.root, recursive = True) if "_build" not in i and "workflow/report" not in i)
        for path in _ipynb:
            # ignore the 'reports' folder name when building the tree
            parts = [i for i in Path(path).parts if i != "reports"]

            if len(parts) == 1:
                # Root level file
                if '__root__' not in tree:
                    tree['__root__'] = []
                tree['__root__'].append(parts[0])
            else:
                # Navigate to parent directory, creating dicts as needed
                current = tree
                for part in parts[:-2]:
                    if part not in current:
                        current[part] = {}
                    current = current[part]

                # Add the file to its parent directory
                parent_dir = parts[-2]

                if parent_dir not in current:
                    current[parent_dir] = set()

                # If it's already a list, append
                if isinstance(current[parent_dir], set):
                    current[parent_dir].add(path)
                # If it's a dict, it has subdirs, so we need to keep it as dict
                # (this handles cases where a dir has both files and subdirs)
        if tree != self.filetree:
            self.filechanges = True
        else:
            self.filechanges = False
        self.filetree = tree

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
    to render them as one navigable and interactive local website, publish it
    to GitHub Pages, etc.
    """

@click.command(context_settings={"help_option_names" : ['-h', '--help']}, epilog = "Documentation: https://pdimens.github.io/harpy/reports/")
@click.option('-g', '--gh-pages', is_flag = True, show_default = True, default = False, help = 'Setup a GitHub Action to build the website on Push')
def init(gh_pages):
    """
    Configure the project for advanced report features
    
    Using MyST, Harpy reports can be built into and rendered as
    an interactive website, which first requires a bit of configuration.
    Use `--gh-pages` to optionally setup a GitHub Action to build
    the website and publish to GitHub on a push to the remote repository.
    """
    notices = []
    if not shutil.which("git"):
        notices.append("- The [green]git[/] software was not found to identify the root directory of this project, so it was assumed this is the root directory.")
        git_dir = ""
    else:
        git_dir = subprocess.run("git rev-parse --show-toplevel".split(), text = True, stdout = subprocess.PIPE).stdout.strip()
    if gh_pages:
        if git_dir:
            fetch_template("buildreports.yml", os.path.join(git_dir, ".github", "workflows", "buildreports.yml"))
        else:
            print_error(
                "not version controlled",
                "Configuring the project and GitHub Action requires this command to be run anywhere within a Git version-controlled directory, however [green]git[/] was unable to detect the root of this repository.",
                "Please verify that this a git-managed repository, and if not, use [blue]git init[/] to set it up as one."
            )
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
@click.option('-r', '--refresh', type = click.IntRange(min = 1, max_open=True), show_default = True, default = 5, help = 'Refresh interval, in seconds')
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
    tracker.overwrite_yaml()
    URL = ""
    try:
        rprint("Starting up the MyST live-server")
        with subprocess.Popen(cmd, cwd = directory, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True) as serve:
            while not URL:
                _myst_output = serve.stdout.readline()
                _url = re.findall(r'http://\S+\s', _myst_output)
                if _url:
                    URL += _url[0].strip()
            rprint(f"MyST live-server started: [blue bold]{URL}[/]\nTerminate it with [bold yellow]ctrl+c[/]")
            while True:
                tracker.scan_for_reports()
                tracker.overwrite_yaml()
                sleep(refresh)
    except KeyboardInterrupt:
        pass

report.add_command(init)
report.add_command(render)



#        def recursive_transform(_toc):
#            """
#            Transform a MyST formatted table of contents into a sensible nested `dict`.
#            """
#            result = {}
#            for item in _toc:
#                title = item['title']
#                children = item['children']
#                
#                # If children is a list of dicts with 'title' key, recurse
#                if isinstance(children, list) and children and isinstance(children[0], dict):
#                    if 'title' in children[0]:
#                        # Nested structure
#                        result[title] = recursive_transform(children)
#                    elif 'file' in children[0]:
#                        # Terminal list of files - extract filenames
#                        result[title] = [child['file'] for child in children]
#                else:
#                    # Children is some other structure (shouldn't happen with our format)
#                    result[title] = children
#            return result