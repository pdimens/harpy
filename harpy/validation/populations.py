from pathlib import Path
from harpy.common.printing import CONSOLE, print_notice, print_error

class Populations():
    '''
    A class to contain and validate a sample-grouping input file.
    '''
    def __init__(self, filename, infiles, quiet:bool = False):
        self.file = filename
        self.quiet: bool = quiet

        if not self.quiet:
            CONSOLE.log("Validating all samples are present between input and populations files")

        rows = []
        with open(self.file, "r", encoding="utf-8") as f:
            for i,j in enumerate(f):
                if j.strip() and not j.lstrip().startswith("#"):
                    rows.append([i,j])
        invalids = [i for i in rows if len(i[-1].split()) < 2]
        if invalids:
            print_error(
                "invalid format",
                f"There are [bold]{len(invalids)}[/] rows in [bold]{self.file}[/] without a space/tab delimiter or don't have two entries for sample[dim]<tab>[/dim]population. Terminating Harpy to avoid downstream errors.",
                f"Make sure every entry in [bold]{self.file}[/] uses space or tab delimeters and has both a sample name and population designation. You may comment out rows with a [green bold]#[/] to have Harpy ignore them.",
                "The rows and values causing this error are",
                [f"{i[0]+1}\t{i[1]}" for i in invalids]
            )

        in_samples = [Path(i).stem for i in infiles]
        popsamples = [i[-1] for i in rows]
        missing_samples = [x for x in popsamples if x not in in_samples]
        overlooked = [x for x in in_samples if x not in popsamples]
        if len(overlooked) > 0:
            print_notice(f"There are [bold]{len(overlooked)}[/] samples found in the inputs that weren\'t included in [blue]{self.file}[/]. This will [bold]not[/] cause errors and can be ignored if it was deliberate. Commenting or removing these lines will avoid this message. The samples are:\n" + ", ".join(overlooked))
        if len(missing_samples) > 0:
            print_error(
                "mismatched inputs",
                f"There are [bold]{len(missing_samples)}[/] samples included in [blue]{self.file}[/] that weren\'t found in in the inputs. Terminating Harpy to avoid downstream errors.",
                f"Make sure the spelling of these samples is identical in the inputs and [blue]{self.file}[/], or remove them from [blue]{self.file}[/].",
                "The samples causing this error are",
                ", ".join(sorted(missing_samples))
            )
