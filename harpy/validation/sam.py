
from itertools import chain
import os
import re
from rich.markdown import Markdown
from harpy.common.printing import print_error

class SAM():
    '''
    A class to contain and validate SAM input files.
    '''
    def __init__(self, filenames):
        if any(isinstance(i, list) for i in filenames):
            self.files = list(chain.from_iterable(filenames))
        else:
            self.files = filenames
        self.count = 0

        re_ext = re.compile(r"\.(bam|sam)$", re.IGNORECASE)
        uniqs = set()
        dupes = []
        inv_pattern = r'[^a-zA-Z0-9._-]+'
        badmatch = []
        for i in self.files:
            bn = os.path.basename(re_ext.sub("", str(i)))
            if bn in uniqs:
                dupes.append(bn)
            else:
                uniqs.add(bn)
                self.count += 1
            if re.search(inv_pattern, os.path.basename(i)):
                badmatch.append(os.path.basename(i))
        
        if badmatch:
            print_error(
                "invalid characters",
                "Invalid characters were detected in the input file names.",
                Markdown("Valid file names may contain only:\n- **A-Z** characters (case insensitive)\n- **.** (period)\n- **_** (underscore)\n- **-** (dash)"),
                "The offending files",
                ", ".join(badmatch)
                )
        if dupes:
            dupe_out = []
            for i in dupes:
                dupe_out.append(" ".join([j for j in self.files if i in j]))
            print_error(
                "clashing sample names",
                Markdown("Identical filenames were detected, which will cause unexpected behavior and results.\n- files with identical names but different-cased extensions are treated as identical\n- files with the same name from different directories are also considered identical"),
                "Make sure all input files have unique names.",
                "Files with clashing names",
                dupe_out
            )

