import json as js
import click
import sys
import pandas as pd

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

def process_json(json):
    pheniqs = js.load(open(json, 'r'))
    d = {}
    total = pheniqs['incoming']['count']
    for i,j in enumerate(["A","B","C"],1):
        #d[f'{j}-unclassified'] = round(1- pheniqs['molecular'][i]['classified fraction'], 5)
        for k in pheniqs['molecular'][i]['classified']:
            barcode = k['barcode'][0]
            if barcode not in d:
                d[barcode] = Molecule()
            d[barcode].update(j, k['pooled fraction'])
    return total,d

@click.command(no_args_is_help = True, epilog = "Documentation: https://pdimens.github.io/harpy/workflows/preprocess/")
@click.argument("json", required = True, type=click.Path(exists = True, dir_okay=False, resolve_path=True))
@click.help_option('--help', hidden = True)
def preproc_post(json):
    total, dd = process_json(json)
    #sys.stdout.write(f"Total: {dd['total']}\n")
    #del dd['total']
    #sys.stdout.write(f"A Unclassified: {dd['A-unclassified']}\n")
    #sys.stdout.write(f"B Unclassified: {dd['B-unclassified']}\n")
    #sys.stdout.write(f"C Unclassified: {dd['C-unclassified']}\n")
    #del dd['A-unclassified']
    #del dd['B-unclassified']
    #del dd['C-unclassified']
    #for k,v in dd.items():
    #    sys.stdout.write(f"{k}\t{v}\n")
    df = pd.DataFrame.from_dict(
        {bc: vars(mol) for bc, mol in dd.items()},
        orient='index'
    ).rename_axis('barcode').reset_index()
    print(df)