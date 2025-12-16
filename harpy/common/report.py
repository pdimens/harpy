"""Harpy module to create a sample grouping file"""

import glob
import os
from os.path import isfile
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
            self.toc_to_filetree(_yml["project"].pop("toc"))
            #_yml["project"].pop("toc")
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
        
    def clean_filetree(self):
        """
        Scans the filetree for nonexistant files or folders without ipynb files
        and removes them, then recusively cleans up keys that terminate in an
        empty value
        """
        def clean_paths(nested_dict, base_path=''):
            """
            Traverse nested dict, build complete paths, and remove entries for 
            non-existent paths or directories with no .ipynb files.

            Returns True if the current dict should be kept, False if it should be removed.
            """
            if '_items' in nested_dict:
                # We're at a terminal node
                items_to_remove = set()

                for item in nested_dict['_items']:
                    # Build complete path
                    if base_path:
                        full_path = os.path.join(base_path, item)
                    else:
                        # For root items, use them as-is
                        full_path = item

                    # Check if path exists and has .ipynb files
                    should_remove = False

                    if not os.path.exists(full_path):
                        should_remove = True
                    elif os.path.isfile(full_path):
                        # For files (like root items), keep them if they exist
                        should_remove = False
                    elif os.path.isdir(full_path):
                        # Check if directory contains any .ipynb files
                        has_ipynb = False
                        for root, dirs, files in os.walk(full_path):
                            if any(f.endswith('.ipynb') for f in files):
                                has_ipynb = True
                                break

                        if not has_ipynb:
                            should_remove = True

                    if should_remove:
                        items_to_remove.add(item)

                # Remove invalid items
                nested_dict['_items'] -= items_to_remove

                # If _items is now empty, this node should be removed
                return len(nested_dict['_items']) > 0

            else:
                # We're at an intermediate node, recurse into children
                keys_to_remove = []

                for key, value in list(nested_dict.items()):
                    # Special handling for 'root' - don't add it to the path
                    if key == 'root':
                        new_base_path = ''
                    else:
                        # Build the path for this level
                        new_base_path = os.path.join(base_path, key) if base_path else key

                    # Recursively clean
                    should_keep = clean_paths(value, new_base_path)

                    if not should_keep:
                        keys_to_remove.append(key)
                # Remove empty branches
                for key in keys_to_remove:
                    del nested_dict[key]
                # Keep this node if it has any remaining children
                return len(nested_dict) > 0

        self.filechanges = clean_paths(self.filetree)

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

        def parse_toc(toc_list):
            result = {}
            
            def add_path_to_dict(path, target_dict):
                """Add a file path to the nested dictionary structure"""
                parts = path.split('/')
                current = target_dict
                
                # Navigate/create nested structure for all but the last part
                for part in parts[:-1]:
                    if part not in current:
                        current[part] = {}
                    current = current[part]
                
                # Add the final directory to _items set
                final_part = parts[-1]
                if '_items' not in current:
                    current['_items'] = set()
                current['_items'].add(final_part)
            
            def process_item(item):
                # If it's a top-level file, add to root
                if 'file' in item:
                    if 'root' not in result:
                        result['root'] = {'_items': set()}
                    result['root']['_items'].add(item['file'])
                    return
                
                # If it has children, process them recursively
                if 'children' in item:
                    for child in item['children']:
                        # If child has a pattern, extract and add the path
                        if 'pattern' in child:
                            pattern = child['pattern']
                            # Remove the glob pattern to get just the directory path
                            path = pattern.split('/*')[0]  # This handles any /** pattern
                            add_path_to_dict(path, result)
                        else:
                            # Continue processing nested items
                            process_item(child)
            
            # Process each top-level item
            for item in toc_list:
                process_item(item)

            return result

        self.filetree = parse_toc(toc)
        self.clean_filetree()
        CONSOLE.print(self.filetree, style = "magenta")
        exit(0)

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
