"""Harpy module to create a sample grouping file"""

import glob
import hashlib
from importlib import resources
import json
import os
from pathlib import Path
import shutil
import subprocess
import uuid
import yaml
from harpy.common.printing import HarpyPrint
from harpy import __version__

class ReportRender():
    def __init__(self, root: str = "", markdown: bool = False):
        """
        Initialize the ReportRender for a given project root and optional Markdown support.
        
        Parameters:
            root (str): Filesystem path to the project root where myst.yml and content are located.
            markdown (bool): When True, include Markdown (.md) files alongside notebooks (.ipynb)
                during scanning and TOC generation.
        
        Behavior:
            - Initializes instance attributes: `root`, `configfile`, empty `filetree`, `print` helper,
              and `markdown`.
            - Loads or creates the MyST configuration by calling `init()`.
            - If the loaded config contains a `project.toc`, converts that TOC into `filetree` via
              `toc_to_filetree()` and removes the `toc` entry from the loaded config.
            - Otherwise computes and stores a checksum for the current `filetree`.
            - Stores the loaded/updated YAML config in `self.config`.
        """
        self.root = root
        self.configfile = os.path.join(root, "myst.yml")
        self.filetree: dict = {}
        self.print = HarpyPrint()
        self.markdown = markdown
        yml = self.init()
        if "toc" in yml["project"]:
            self.toc_to_filetree(yml["project"].pop("toc"))
        else:
            self.checksum = checksum(self.filetree)
        self.config: dict = yml

    def init(self):
        """
        Initialize and prepare the MyST working directory and configuration.
        
        Ensures the project's myst.yml exists (creates and writes a default when missing) and loads it. Creates a .report directory with an index.md (templated with the package version) and ensures favicon.ico, logo.svg, and logo-dark.svg exist by copying packaged assets when absent. Validates that the loaded YAML is a mapping containing top-level "project" and "site" keys and emits an error via self.print.error if validation fails.
        
        Returns:
            dict: The loaded or newly created MyST configuration mapping.
        """
        if not os.path.exists(self.configfile):
            _yml = myst_yaml()
            self.scan()
            with open(self.configfile, "w") as yml:
                yaml.dump(_yml, yml, default_flow_style= False, sort_keys=False, width=float('inf'))
        else:
            with open(self.configfile, "r") as yml:
                _yml = yaml.full_load(yml)

        if not os.path.isfile(os.path.join(self.root, ".report", "index.md")):
            os.makedirs(os.path.join(self.root, ".report"), exist_ok=True)
            with (
                resources.as_file(resources.files("harpy.templates") / "report_index.md") as _source,
                open(_source, 'r') as md,
                open(os.path.join(self.root, ".report", "index.md"), 'w') as indexmd,
                ):
                indexmd.write(md.read().format(__version__))

        if not os.path.isfile(os.path.join(self.root, ".report", "favicon.ico")):
            shutil.copy(str(resources.files("harpy.report") / 'favicon.ico'), os.path.join(self.root, ".report", "favicon.ico"))
        if not os.path.isfile(os.path.join(self.root, ".report", "logo.svg")):
            shutil.copy(str(resources.files("harpy.report") / 'logo.svg'), os.path.join(self.root, ".report", "logo.svg"))
        if not os.path.isfile(os.path.join(self.root, ".report", "logo-dark.svg")):
            shutil.copy(str(resources.files("harpy.report") / 'logo-dark.svg'), os.path.join(self.root, ".report", "logo-dark.svg"))

        if not isinstance(_yml, dict) or "project" not in _yml or "site" not in _yml:
            self.print.error(
                "invalid MyST configuration",
                "The [blue]myst.yml[/] file that was found is improperly formatted.",
                "If you manually edited this file, please try to remake it using [blue]harpy report init[/], otherwise check that the file is in proper YAML format and has keys for [green]project[/] and [green]site[/]."
            )
        return _yml

    def update_yaml(self):
        """
        Update the MyST config file with a TOC generated from the current filetree if the filetree changed.
        
        If the computed checksum of the current filetree differs from the stored checksum, this method generates a MyST-compatible table of contents from self.filetree, writes the updated YAML to self.configfile, removes the temporary `"toc"` key from self.config["project"], and updates self.checksum to the new value.
        """
        if self.checksum == checksum(self.filetree):
            return
        self.config["project"]["toc"] = self.filetree_to_toc()
        with open(self.configfile, "w") as yml:
            yaml.dump(self.config, yml, default_flow_style= False, sort_keys=False, width=float('inf'))
        del self.config["project"]["toc"]
        self.checksum = checksum(self.filetree)

    def scan(self):
        """
        Builds and stores a nested file-tree of project content found under self.root.
        
        Searches recursively for files ending with `.ipynb` and, when `self.markdown` is True, `.md`, excluding any paths that contain `_build` or start with `workflow/`. Populates `self.filetree` as a nested dictionary representing directory hierarchy, where each directory node may contain a `_items` set of terminal names found in that directory.
        """
        _ipynb = set(i for i in glob.iglob("**/*.ipynb", root_dir = self.root, recursive = True) if "_build" not in i and "workflow/" not in i)
        if self.markdown:
            _md = set(i for i in glob.iglob("**/*.md", root_dir = self.root, recursive = True) if "_build" not in i and "workflow/" not in i)
            _ipynb = _ipynb.union(_md)
        # filter out ""
        _dirs = set(d for d in (os.path.dirname(i) for i in _ipynb) if d)
        for path in _dirs:
            parts = Path(path).parts
            current = self.filetree           
            # Navigate/create nested structure
            for part in parts[:-1]:
                if part not in current:
                    current[part] = {}
                current = current[part]

            # Add final directory to a set
            final_dir = parts[-1]
            if '_items' not in current:
                current['_items'] = set()
            current['_items'].add(final_dir)

    def clean_filetree(self):
        """
        Prune self.filetree by removing entries that no longer exist or point to directories without .ipynb or .md files.
        
        Traverses the nested filetree in-place, removes terminal items whose filesystem path is missing or is a directory that contains no .ipynb or .md files (case-insensitive), and then recursively deletes intermediate dictionary nodes that become empty as a result. Emits a dimmed removal message for each removed path via self.print.
        """
        def clean_paths(nested_dict, base_path=''):
            """
            Prune a nested filetree node in-place by removing entries whose paths do not exist or directories that contain no .ipynb or .md files.
            
            This function operates recursively on a nested_dict that either contains a terminal '_items' set of file/directory names or child dictionaries representing subdirectories. It removes invalid entries from '_items' or deletes empty child keys, emits a removal message via self.print.print for each removed path, and mutates nested_dict in-place.
            
            Parameters:
            	nested_dict (dict): A node in the filetree structure; terminal nodes contain a set at key '_items'.
            	base_path (str): The path prefix to prepend when resolving items in this node (empty for root).
            
            Returns:
            	True if the node contains remaining items or children and should be kept, `False` otherwise.
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
                            if any(f.lower().endswith('.ipynb') or f.lower().endswith('.md') for f in files):
                                has_ipynb = True
                                break

                        if not has_ipynb:
                            should_remove = True

                    if should_remove:
                        self.print.print(f"[removed - no md/ipynb files] " + os.path.relpath(full_path), style = "dim", markup = False, soft_wrap = True)
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
        clean_paths(self.filetree)

    def filetree_to_toc(self):
        """
        Convert the ReportRender.filetree into a MyST-compatible table of contents structure.
        
        The returned list contains TOC entries as dictionaries. Possible entry shapes:
        - {"file": "<filename>"} for individual files.
        - {"title": "<name>", "children": [...]} for nested sections.
        - {"pattern": "<glob>"} for directory globs that match notebook files (and, when markdown is enabled, markdown files).
        
        Returns:
            list: A MyST `toc`-formatted list of dictionaries representing files, sections, and glob patterns.
        """
        def recursive_transform(tree, origpath = ""):
            """
            Convert a nested filetree structure into a MyST-compatible table-of-contents (TOC) list.
            
            Parameters:
                tree (dict | set | list): Nested filetree node(s) produced by ReportRender.scan()/toc parsing. Nodes may contain:
                    - dicts representing directories with optional "root" -> {"_items": [...]}
                    - sets of item names
                    - non-dict leaf values treated as directory leaves.
                origpath (str): Base path used to build pattern globs for leaf directories (defaults to "").
            
            Returns:
                list[dict]: A list of TOC entries where each entry is one of:
                    - {'file': <filename>} for direct file entries,
                    - {'title': <title>, 'children': [...]} for nested sections,
                    - {'pattern': <glob>} for directory globs (e.g., "<path>/**.ipynb" and, if markdown enabled, "<path>/**.md").
            """
            build_path = origpath
            if isinstance(tree, dict):
                result = []
                [result.append({"file" : i}) for i in tree.get("root", {}).get("_items", [])]

                for key, value in sorted(tree.items(), key=lambda x: x[0].casefold()):
                    if key == "root":
                        continue
                    if isinstance(value, dict):
                        # Recursively transform the nested dict, but skip over "reports"
                        if key == "reports":
                            return recursive_transform(value, os.path.join(build_path, key))
                        result.append({
                            'title': key,
                            'children': recursive_transform(value, os.path.join(build_path, key))
                        })
                    elif isinstance(value, set):
                        for i in value:
                            if "reports" in i:
                                result.append({'pattern': f"{os.path.join(origpath, i, '**.ipynb')}"})
                                if self.markdown:
                                    result.append({'pattern': f"{os.path.join(origpath, i, '**.md')}"})
                            else:
                                result.append({
                                    'title': i,
                                    'children': recursive_transform(value, os.path.join(origpath, i))
                                })
                return result
            else:
                # Value is a list (leaf node) - convert directory name to {'pattern': '*.ipynb}
                wildcard = os.path.join(origpath, "**.ipynb")
                if self.markdown:
                    return [{'pattern': wildcard}, {'pattern': os.path.join(origpath, "**.md")}]
                return [{'pattern': wildcard}]

        return recursive_transform(self.filetree)

    def toc_to_filetree(self, toc):
        """
        Parse a MyST TOC list into the renderer's nested filetree.
        
        Converts a TOC (list of mapping items as used in MyST/myst.yml) into the internal
        nested dictionary structure stored at self.filetree, computes and stores a checksum
        for that structure, and prunes any paths that do not exist on disk.
        
        Parameters:
            toc (list): A MyST-compatible TOC represented as a list of mapping items
                (each item may contain keys like "file", "children", or "pattern").
        """
        def parse_toc(toc_list):
            """
            Parse a MyST TOC list into a nested filetree dictionary used by ReportRender.
            
            Parameters:
                toc_list (list): A MyST `toc` list where items may contain `file`, `children`, and `pattern` entries.
            
            Returns:
                dict: Nested dictionary representing the project layout. Directory keys map to sub-dictionaries; terminal entries are stored in sets under the `_items` key. Top-level `file` entries are added to `root`->_items. `pattern` entries are interpreted as directory paths (the portion before any `/*` glob).
            """
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
        self.checksum = checksum(self.filetree)
        self.clean_filetree()        

class StableEncoder(json.JSONEncoder):
    def default(self, o):
        """
        Provide JSON-serializable representations for objects not handled by the base encoder.
        
        When `o` is a `set`, return a dict `{"__set__": [...]}` where the list contains the set elements in sorted order to ensure deterministic serialization. For any other object, defer to the base JSONEncoder's `default` implementation.
        
        Parameters:
            o: The object to encode.
        
        Returns:
            A JSON-serializable representation of `o`; for sets, a dict with key `__set__` and a sorted list of elements.
        """
        if isinstance(o, set):
            return {"__set__": sorted(o)}  # sorted() ensures stable ordering
        return super().default(o)

# Source - https://stackoverflow.com/a/66583613
# Posted by Chris Maes, modified to serialize sets
# Retrieved 2026-05-19, License - CC BY-SA 4.0
def checksum(d):
    """
    Compute an MD5 checksum of a Python object by serializing it to a stable JSON string.
    
    Uses StableEncoder to deterministically serialize `set` objects, sorts object keys, and ensures ASCII-safe encoding before computing the MD5 digest.
    
    Parameters:
        d (Any): The Python object (typically a dict) to serialize and checksum.
    
    Returns:
        str: Hexadecimal MD5 digest of the serialized JSON representation.
    """
    return hashlib.md5(
        json.dumps(d, cls=StableEncoder, sort_keys=True, ensure_ascii=True).encode('utf-8')
    ).hexdigest()

def rand_id() -> str:
    """
    Generate a hyphenated random UUID string formatted as 8-4-4-4-12.
    
    Returns:
        str: UUID string composed of hexadecimal characters in five groups (8-4-4-4-12), 36 characters including hyphens.
    """
    uid = uuid.uuid4().hex
    indices = [0, 8, 12, 16, 20, 32]
    return "-".join(uid[i:j] for i,j in zip(indices, indices[1:]))

def myst_yaml() -> dict:
    """
    Constructs a default MyST configuration dictionary for a Harpy report site.
    
    Attempts to read the repository's origin URL and, if available, includes it as the project's `github` field. The returned dictionary contains top-level keys `version`, `site` (theme, actions, and UI options such as logo and favicon), and `project` (an `id`, optional `github`, `edit_url`, `title`, `description`, and an initial `toc` pointing to `.report/index.md`).
    
    Returns:
        dict: A MyST-compatible configuration dictionary for use as a default `myst.yml`.
    """
    try:
        git_url = subprocess.run("git remote get-url origin".split(), text = True, stdout = subprocess.PIPE).stdout.strip().removesuffix(".git")
    except Exception:
        git_url = ""
    return {
        "version" : 1,
        "site" : {
            "template": "book-theme",
            "actions" : [
                {"title": "📖 Docs", "url": "https://pdimens.github.io/harpy"}
            ],
            "options" : {
                "favicon" : ".report/favicon.ico",
                "folders": True,
                "hide_footer_links" : True,
                #"hide_myst_branding" : True,
                "logo" : ".report/logo.svg",
                "logo_dark" : ".report/logo-dark.svg",
                "logo_text" : "Harpy Reports"
            },
        },
        "project" : {
            "id": rand_id(),
            **({"github" : git_url} if git_url else {}),
            "edit_url": 'null',
            "title" : f"{os.path.basename(os.getcwd())}",
            "description" : "Harpy reports, powered by Jupyter and MyST.",
            "toc": [{"file" : ".report/index.md"}]
        }
    }
