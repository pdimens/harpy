import os
import sys
import rich_click as click
import shutil
import subprocess
import uuid
from harpy.common.file_ops import fetch_template
from harpy.common.printing import print_error, print_notice, CONSOLE

def rand_id():
    uid = uuid.uuid4().hex
    indices = [0, 8, 12, 16, 20, 32]
    return "-".join(uid[i:j] for i,j in zip(indices, indices[1:]))

def myst_yaml() -> str:
    git_url = subprocess.run("git remote get-url origin".split(), text = True, stdout = subprocess.PIPE).stdout.strip().removesuffix(".git")
    gh_line = f"\n  github: {git_url}" if git_url else ""
    return f"""# See docs at: https://mystmd.org/guide/frontmatter
version: 1
site:
  template: book-theme
  # options:
  #   favicon: favicon.ico
  #   logo: site_logo.png
project:
  id: {rand_id()}
  title: Harpy Reports
  description: The various reports produced by Harpy, aggregated into a navigable website using MyST.
  authors: ["Pavel Dimens"]{gh_line}
"""

@click.group(context_settings={"help_option_names" : []})
def reports():
    """
    Setup and render harpy reports

    Harpy reports are provided as Jupyter Notebooks. The subsequent subcommands
    will configure your project directory to make these much nicer, enabling you
    to render them as one navigable and interactive local website, publish it
    to GitHub, etc.
    """

@click.command(epilog = "Documentation: https://pdimens.github.io/harpy/reports/")
@click.option('-g', '--gh-pages', is_flag = True, show_default = True, default = False, help = 'Setup a GitHub Action to build the website on Push')
def init(gh_pages):
    """
    Configure the repository for advanced report features
    
    Using MyST, the harpy reports can be rendered and built into
    an interactive website, which first requires a bit of configuration.
    Use `--gh-pages` to optionally setup a GitHub Action to build
    the website and publish to GitHub on a push to the remote repository.
    """
    notices = []
    if not shutil.which("git"):
        notices.append("- The [green]git[/] software was not found to identify the root directory of this project, so it was assumed this is the root directory.")
        git_dir = ""
    else:
        git_dir = subprocess.run("git rev-parse --show-toplevel".split(), text = True, stdout = subprocess.PIPE).stdout.strip()
    if gh_pages:
        if git_dir:
            fetch_template("buildreports.yml", os.path.join(git_dir, ".github", "workflows", "buildreports.yml"))
        else:
            print_error(
                "not version controlled",
                "Configuring the project and GitHub Action requires this command to be run anywhere within a Git version-controlled directory, however [green]git[/] was unable to detect the root of this repository.",
                "Please verify that this a git-managed repository, and if not, use [blue]git init[/] to set it up as one."
            )
    #TODO pull any specific configurations here too
    with open(os.path.join(git_dir, "myst.yml"), "w") as yml:
        yml.write(myst_yaml())
    if not shutil.which("myst"):
        notices.append("- The [green]MyST[/] software is required for a locally-served report website, however it was not found in this environment.")
    if notices:
        print_notice("\n".join(notices))
reports.add_command(init)