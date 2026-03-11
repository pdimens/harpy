import click
import sys
import pandas as pd
import pysam
import re

class Mean():
    def __init__(self):
        self.sum = 0.0
        self.n = 0
    def add(self, val):
        self.sum += val
        self.n += 1

class Molecule():
    def __init__(self):
        self.A = 0.0
        self.B = 0.0
        self.C = 0.0
    def update(self, segment, frac: float):
        setattr(self, segment, round(frac,5))
    def __str__(self) -> str:
        return f"{self.A}\t{self.B}\t{self.C}"

reBX = re.compile(r'BX:Z:(\S+)')
reOX = re.compile(r'OX:Z:(\S+)')

@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/preprocess/")
@click.argument("input", required = True, type=click.Path(exists = True, dir_okay=False, resolve_path=True))
@click.help_option('--help', hidden = True)
def preproc_post(input):
    with pysam.FastxFile(input, persist=False) as fq:
        for record in fq:
            comments = record.comment.split()
            BX 
