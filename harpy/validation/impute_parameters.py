
import re
from rich.table import Table
from rich import box
from harpy.common.printing import print_error

class ImputeParams():
    '''
    A class to contain and validate a STITCH imputation parameter file.
    Validation checks the STITCH parameter file for column names, order, types, missing values, etc.
    '''
    def __init__(self, filename):
        self.file: str = filename

        with open(self.file, "r", encoding="utf-8") as paramfile:
            header = paramfile.readline().rstrip().lower()
            headersplt = header.split()
            colnames = ["name", "model", "usebx", "bxlimit", "k", "s", "ngen"]
            correct_header = sorted(colnames)
            n_cols = len(colnames)
            row = 1
            badrows = []
            badlens = []
            if sorted(headersplt) != correct_header:
                invalid_col = [i for i in headersplt if i not in colnames]
                missing_col = [i for i in colnames if i not in headersplt]
                errtext = []
                culprit_text = []
                if invalid_col:
                    errtext.append("\n  - invalid column names")
                    culprit_text.append("[red]Invalid columns:[/] " + ", ".join(invalid_col))
                if missing_col:
                    errtext.append("\n  - missing columns")
                    culprit_text.append("[yellow]Missing columns:[/] " + ", ".join(missing_col))

                print_error(
                    "incorrect columns",
                    f"Parameter file [bold]{self.file}[/] has incorrect column names" + "".join(errtext) + "\nValid names are: [green bold]" + " ".join(colnames) + "[/]",
                    f"Fix the headers in [bold]{self.file}[/] or use [blue bold]harpy template[/] to generate a valid parameter file and modify it with appropriate values.",
                    "Column causing this error",
                    "\n".join(culprit_text)
                )
            else:
                # in case the columns are out of order, reorder the columns to match `colnames`
                col_order = [colnames.index(item) for item in headersplt]

            # instantiate dict with colnames
            self.parameters: dict = {}
            self.count: int = 0
            while True:
                # Get next line from file
                line = paramfile.readline()
                row += 1
                # if line is empty, end of file is reached
                if not line:
                    break
                if line == "\n":
                    # be lenient with empty rows
                    continue
                # split the line by whitespace and reorder to match expected colname order
                row_values = [line.rstrip().split()[i] for i in col_order]
                if len(row_values) == n_cols: 
                    self.parameters[row_values[0]] = dict(zip(colnames[1:], row_values[1:]))
                    self.count += 1
                else:
                    badrows.append(row)
                    badlens.append(len(row_values))
            if len(badrows) > 0:
                _outrows = []
                for i in zip(badrows, badlens):
                    _outrows.append(f"{i[0]}\t{i[1]}")
                print_error(
                    "invalid rows",
                    f"Parameter file [blue]{self.file}[/] is formatted incorrectly. Not all rows have the expected {n_cols} columns.",
                    f"See the problematic rows below. Check that you are using a whitespace (space or tab) delimeter in [blue]{self.file}[/] or use [blue green]harpy template[/blue green] to generate a valid parameter file and modify it with appropriate values.",
                    "Rows causing this error and their column count",
                    _outrows
                )
            
            # validate each row
            row_error = False
            errtable = Table(show_footer=True, box=box.SIMPLE)
            errtable.add_column("Row", justify="right", style="white", no_wrap=True)
            errtable.add_column("Columns with Issues", style = "white")

            row = 1
            for k,v in self.parameters.items():
                badcols = []
                if re.search(r'[^a-z0-9\-_\.]',k, flags = re.IGNORECASE):
                    badcols.append("name")
                if v["model"] not in ["pseudoHaploid", "diploid","diploid-inbred"]:
                    badcols.append("model")
                if v["usebx"].lower() not in ["true", "false", "yes", "y", "no", "n"]:
                    badcols.append("usebx")
                else:
                    if v["usebx"].lower() in ["true", "yes", "y"]:                
                        v["usebx"] = True
                    else:
                        v["usebx"] = False
                for param in ["bxlimit", "k", "s", "ngen"]:
                    if not v[param].isdigit():
                        badcols.append(param)
                    else:
                        v[param] = int(v[param])
                if badcols:
                    row_error = True
                    errtable.add_row(f"{row}", ", ".join(badcols))
                row += 1
            if row_error:
                print_error(
                    "invalid parameter values",
                    f"Parameter file [bold]{self.file}[/] is formatted incorrectly. Some rows have incorrect values for one or more parameters.",
                    "Review the table below of which rows/columns are causing issues",
                    "Formatting Errors",
                    errtable
                )