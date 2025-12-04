"""Harpy module to create a sample grouping file"""

import copy
import glob
import os
from pathlib import Path
import uuid
import shutil
import subprocess
import yaml
from harpy.common.file_ops import fetch_template
from harpy.common.printing import print_error, print_notice

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
