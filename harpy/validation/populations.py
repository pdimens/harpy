from pathlib import Path
from harpy.common.printing import print_notice, print_error, print_solution_offenders

class Populations():
    '''
    A class to contain and validate a sample-grouping input file.
    '''
    def __init__(self, filename, infiles):
        self.file = filename
        
        with open(self.file, "r", encoding="utf-8") as f:
            popsamples = [i.split()[0] for i in f.readlines() if i != "\n" and not i.lstrip().startswith("#")]
        in_samples = [Path(i).stem for i in infiles]
        missing_samples = [x for x in popsamples if x not in in_samples]
        overlooked = [x for x in in_samples if x not in popsamples]
        if len(overlooked) > 0:
            print_notice(f"There are [bold]{len(overlooked)}[/] samples found in the inputs that weren\'t included in [blue]{self.file}[/]. This will [bold]not[/] cause errors and can be ignored if it was deliberate. Commenting or removing these lines will avoid this message. The samples are:\n" + ", ".join(overlooked))
        if len(missing_samples) > 0:
            print_error(
                "mismatched inputs",
                f"There are [bold]{len(missing_samples)}[/] samples included in [blue]{self.file}[/] that weren\'t found in in the inputs. Terminating Harpy to avoid downstream errors.",
                False
            )
            print_solution_offenders(
                f"Make sure the spelling of these samples is identical in the inputs and [blue]{self.file}[/], or remove them from [blue]{self.file}[/].",
                "The samples causing this error are",
                ", ".join(sorted(missing_samples))
            )
