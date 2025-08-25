"""
Various system-level checks such as whether packages exist, architecture, etc.
"""

import os
import platform
import shutil
from harpy.common.printing import print_error, print_notice
import rich_click as click
from rich.markdown import Markdown

def is_conda_package_installed(package_name) -> int:
    """
    Return 0 if package exists
    Return 1 if conda is not found but CONDA_PREFIX is defined, implying it's mamba
    Return 2 if package not found in current conda env
    Return 3 if error
    Return 4 if CONDA env not detected
    """
    if "CONDA_PREFIX" in os.environ:
        try:
            from conda.core.prefix_data import PrefixData
            from conda.base.context import context
        except ModuleNotFoundError:
            # CONDA_PREFIX is there but conda itself isnt: likely mamba
            return 1
        try:
            # Get the current environment prefix
            prefix = context.active_prefix or context.root_prefix

            # Get prefix data for current environment
            prefix_data = PrefixData(prefix)

            # Check if package exists in the prefix
            if any(package_name in record.name for record in prefix_data.iter_records()):
                return 0
            else:
                return 2
        except Exception:
            return 3
    else:
        return 4

def is_pip_package_installed(package_name) -> bool:
    from importlib.metadata import PackageNotFoundError
    import importlib.metadata
    try:
        importlib.metadata.version(package_name)
        return True
    except PackageNotFoundError:
        return False

def is_pixi_shell() -> bool:
    # Check for pixi-specific environment variables
    pixi_indicators = [
        'PIXI_PROJECT_ROOT',
        'PIXI_PROJECT_NAME',
        'PIXI_PROJECT_MANIFEST',
        'PIXI_ENVIRONMENT_NAME'
    ]

    return any(var in os.environ for var in pixi_indicators)

def package_absent(pkg: str, executor: bool = True) -> bool:
    """helper function to search for a package in the active conda environment"""
    if executor:
        out_text = "Using this scheduler requires installing a Snakemake plugin which wasn't detected in this environment. "
    else:
        out_text = "Using the `--container` option requires `apptainer`, which wasn't detected in this environment. "

    out_text += "It can be installed with:"

    # check for conda/mamba
    if shutil.which("conda") or shutil.which("mamba"):
        conda_check = is_conda_package_installed(pkg)
        if conda_check == 0:
            return False
        elif conda_check == 2:
            if is_pixi_shell():
                out_text += f"\n\n```bash\npixi add {pkg}\n```"
            else:
                out_text += f"\n\n```bash\nconda install -c bioconda {pkg}\n```"
        elif conda_check == 1:
            out_text += f"\n\n```bash\nmamba install -c bioconda {pkg}\n```"
        if conda_check in [1,2]:
            if executor:
                print_notice(Markdown(out_text))
            else:
                print_error("missing required package", Markdown(out_text))
            return True

    if not is_pip_package_installed(pkg):
        out_text += f"\n\n```bash\npip install {pkg}\n```"
        if executor:
            print_notice(Markdown(out_text))
        else:
            print_error("missing required package", Markdown(out_text))
        return True
    
    return False

def container_ok(ctx, param, value) -> bool:
    """ 
    Check if the system is linux or has apptainer installed
    """
    if value:
        if os.sys.platform != 'linux':
            raise click.BadParameter(
                "Snakemake uses Apptainer (formerly Singularity) to manage containers, which is only available for Linux systems.", ctx, param
            )
        if shutil.which("apptainer"):
            return value
        else:
            raise click.BadParameter(
                "Container software management requires apptainer, which wasn't detected in this environment.", ctx, param
            )
    return value

def is_arm(allowed: bool) -> None:
    """
    Check if the system uses ARM architecture, for which some workflows or their deps might not work. 
    """
    arch = platform.machine().lower()
    # Common ARM identifiers across Linux/macOS: armv7l, armv8l, arm64, aarch64
    _arm = arch.startswith(("arm", "aarch64")) or "arm64" in arch
    if not _arm or allowed:
        return
    else:
        print_error(
            "incompatible architecture",
            "Your system was detected to use an ARM-based processor, for which not all of the software required by this workflow is compatible with and therefore conda will be unable to install. If possible, please try this command again on an x86-based machine (i.e. [blue]Intel[/] or [red]AMD[/] processors)."
        )