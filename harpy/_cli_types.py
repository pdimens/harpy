"""Module with helper function to set up Harpy workflows"""

import os
import click

class IntList(click.ParamType):
    """A class for a click type which accepts an arbitrary number of integers, separated by a comma."""
    name = "int_list"
    def __init__(self, max_entries):
        super().__init__()
        self.max_entries = max_entries

    def convert(self, value, param, ctx):
        try:
            parts = [i.strip() for i in value.split(',')]
            if len(parts) != self.max_entries:
                raise ValueError
            for i in parts:
                try:
                    int(i)
                except:
                    raise ValueError
            return [int(i) for i in parts]
        except ValueError:
            self.fail(f"{value} is not a valid list of integers. The value should be {self.max_entries} integers separated by a comma.", param, ctx)

class KParam(click.ParamType):
    """A class for a click type which accepts any number of odd integers separated by a comma, or the word auto."""
    name = "k_param"
    def convert(self, value, param, ctx):
        try:
            if value == "auto":
                return value
            parts = [i.strip() for i in value.split(',')]
            for i in parts:
                if int(i) % 2 == 0 or int(i) > 128:
                    raise ValueError
            return [int(i) for i in parts]
        except ValueError:
            self.fail(f"{value} is not 'auto' or odd integers <128 separated by a comma.", param, ctx)

class ContigList(click.ParamType):
    """A class for a click type which accepts a file of contigs or a list of contigs separated by a comma."""
    name = "contig_list"
    def convert(self, value, param, ctx):
        # check if it's a file
        if os.path.exists(value):
            if not os.path.isfile(value):
                self.fail(f"{value} is not a file.", param, ctx)
            with open(value, "r") as cont_in:
                return [i.strip() for i in cont_in.readlines()]
        else:
            return [i.strip() for i in value.split(',')]

