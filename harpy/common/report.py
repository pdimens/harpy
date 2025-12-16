"""Harpy module to create a sample grouping file"""

import glob
import os
from pathlib import Path
import uuid
import subprocess
import yaml
from harpy.common.file_ops import fetch_template
from harpy.common.printing import CONSOLE, print_error

class ReportRender():
    def __init__(self, root: str = ""):
        self.root = root
        self.configfile = os.path.join(root, "myst.yml")
        self.filetree: dict = {}
        self.filechanges: bool = False

        if not os.path.exists(self.configfile):
            _yml = myst_yaml()
            self.scan_for_reports()
            print(self.filetree)
            with open(self.configfile, "w") as yml:
                yaml.dump(_yml, yml, default_flow_style= False, sort_keys=False, width=float('inf'))
        else:
            with open(self.configfile, "r") as yml:
                _yml = yaml.full_load(yml)

        if not os.path.isfile(os.path.join(root, ".report", "index.md")):
            fetch_template("report_index.md", os.path.join(root, ".report", "index.md"))
        if not os.path.isfile(os.path.join(root, ".report", "favicon.svg")):
            fetch_template("favicon.png", os.path.join(root, ".report", "favicon.png"))

        if not isinstance(_yml, dict) or "project" not in _yml:
            print_error(
                "invalid MyST configuration",
                "The [blue]myst.yml[/] file that was found is improperly formatted.",
                "If you manually edited this file, please try to remake it using [blue]harpy report init[/], otherwise check that the file is in proper YAML format and has a [green]project[/] section."
            )
        if "toc" in _yml["project"]:
            #self.toc_to_filetree(_yml["project"].pop("toc"))
            _yml["project"].pop("toc")
        self.config: dict = _yml

    def update_yaml(self):
        """
        Convert self.filelist to a MyST-formatted table of contents and overwrite `self.configfile` with the updated table of contents
        """
        import sys
        if not self.filechanges:
            return
        self.config["project"]["toc"] = self.filetree_to_toc()
        #with open(sys.stdout, "w") as yml:
        #with open(self.configfile, "w") as yml:
        #    yaml.dump(self.config, yml, default_flow_style= False, sort_keys=False, width=float('inf'))
        yaml.dump(self.config, sys.stdout, default_flow_style= False, sort_keys=False, width=float('inf'))
        del self.config["project"]["toc"]
        self.filechanges = False

    def scan_for_reports(self):
        """
        Recursively search `self.root` for files ending in `.ipynb`, filtering out those found in the _build/ directory
        and converts that list of file paths into a nested dictionary tree structure stored as `self.filetree`. 
        """
        _ipynb = set(i for i in glob.iglob("**/*.ipynb", root_dir = self.root, recursive = True) if "_build" not in i and "workflow/report" not in i)
        _dirs = set(os.path.dirname(i) for i in _ipynb)
        result = self.filetree
        
        for path in _dirs:
            parts = Path(path).parts
            current = result           
            # Navigate/create nested structure
            for part in parts[:-1]:
                if part not in current:
                    current[part] = {}
                current = current[part]
            
            # Add final directory to a set
            final_dir = parts[-1]
            if '_items' not in current:
                current['_items'] = set()
            if path not in current['_items']:
                self.filechanges = True
            current['_items'].add(final_dir)
        
        CONSOLE.print("FILETREE AFTER SCAN", self.filetree)
        #exit(0)

    def clean_filetree(self):
        """
        Scans the filetree for nonexistant files and removes them, then recusively cleans up keys that terminate
        in an empty value
        """
        def remove_nonexistent_files(d):
            cleaned = {}
            for key, value in d.items():
                if isinstance(value, dict):
                    # Recursively process nested dictionaries
                    cleaned_files = remove_nonexistent_files(value)
                    # Only keep if the cleaned dict is not empty
                    if cleaned_files:
                        cleaned[key] = cleaned_files
                else:
                    if isinstance(value, set):
                        _passed = set()
                        for _dir in value:
                            empty = True
                            for _ in os.scandir(_dir):
                                _passed.add(_dir)
                                break
                            if _passed:
                                cleaned[key] = _passed
                    elif value not in (None, '', [], {}, (), set()):
                        cleaned[key] = value 
            return cleaned

        _ = remove_nonexistent_files(self.filetree)
        if self.filetree != _:
            self.filechanges = True
            self.filetree = _

    def filetree_to_toc(self):
        """
        Recursively parses `self.filetree` to format it as a MyST table of contents and return it
        """
        def recursive_transform(tree, origpath = ""):
            """
            Transform a nested dict into a format where dict values become 
            {'title': key, 'children': transformed_value}.
            Terminal lists of filenames become lists of {'file': filename} dicts.
            
            Returns a transformed structure with title/children format
            """
            build_path = origpath
            if isinstance(tree, dict):
                result = []

                for key, value in sorted(tree.items()):
                    #build_path = os.path.join(build_path)
                    if isinstance(value, dict):
                        # Recursively transform the nested dict
                        result.append({
                            'title': key,
                            'children': recursive_transform(value, os.path.join(build_path, key))
                        })
                    elif isinstance(value, set):
                        # Value is a terminal directory
                        for i in value:
                            if "reports" in i:
                                result.append({'pattern' : f"{os.path.join(origpath, "reports", "**.ipynb")}"})
                            else:
                                result.append({
                                    'title': i,
                                    'children': recursive_transform(value, os.path.join(origpath, i))
                                })
                    elif isinstance(value, str):
                        result.append({
                            'title': value,
                            'potato': os.path.join(origpath, key)
                        })
                return result
            else:
                # Value is a list (leaf node) - convert directory name to {'pattern': '*.ipynb}
                wildcard = os.path.join(origpath, "**.ipynb")
                return [{'pattern': wildcard}]


        #yaml.dump(recursive_transform(self.filetree))
        return recursive_transform(self.filetree)
        #CONSOLE.print("\n[yellow]TOC[/]", recursive_transform(self.filetree))
        #exit(0)

    def toc_to_filetree(self, toc):
        """
        Recursively parses a `toc` (as interpreted in myst.yml) to format it as a nested `dict` stored in `self.filetree`.
        Removes path trees that don't exist.
        """
        def recursive_transform(transformed):
            if isinstance(transformed, list):
                result = {}
                root_files = set()

                for item in transformed:
                    # Handle case where item is just {"file": filename} at top level
                    if 'file' in item and 'title' not in item:
                        root_files.add(item['file'])
                        continue

                    title = item['title']
                    children = item.get('children', None)
                    if not children:
                        children = item.get('pattern', None)

                    if not children:
                        continue

                    # If children is a list of dicts with 'title' key, make recursive
                    if isinstance(children, list) and isinstance(children[0], dict):
                        if 'title' in children[0]:
                            # Nested structure
                            result[title] = recursive_transform(children)
                        elif 'file' in children[0]:
                            # Terminal list of files - extract filenames
                            result[title] = set(child['file'] for child in children)
                        elif isinstance(children, str):
                            if title not in result:
                                result[title] = set()
                            result[title].add(children)
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
        self.clean_filetree()

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
    try:
        git_url = subprocess.run("git remote get-url origin".split(), text = True, stdout = subprocess.PIPE).stdout.strip().removesuffix(".git")
    except Exception:
        git_url = ""
    d = {
        "version" : 1,
        "site" : {
            "template": "book-theme",
            "actions" : [{"title": "Documentation", "url": "https://pdimens.github.io/harpy"}],
            "options" : {
                "favicon" : ".report/favicon.png",
                "logo" : ".report/favicon.png",
                "logo_text" : "Harpy Reports",
                "hide_footer_links" : True,
                "hide_myst_branding" : True,
                "folders": True
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
