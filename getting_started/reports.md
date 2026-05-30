# :icon-graph: Harpy Reports

Harpy has always maintained a robust reporting system for most workflows. These reports detail the common
characteristics of data after a given workflow (e.g. alignment or phasing), including values, tables, and
figures. The reports were originally written in R/RMarkdown, then ported to Quarto. Once the limits of Quarto
were reached, Harpy 4.0 introduced a complete overhaul of the reporting system using Jupyter Notebooks. With
Jupyter comes several benefits:
- Code and output are stored in the notebook
- GitHub, JupyterLab, and VScode (and derivatives) natively render notebooks nicely
- Harpy can leverage [MyST](https://mystmd.org/) (via Jupyter Book) to render everything into a _cohesive_ report webiste

!!!
Complete overhaul is not an overstatement-- all the R code was ported into Python and reformatted for Jupyter format.
That meant 100% new code **everywhere**, all new plots, completely new machinery for local rendering, an entirely
new `harpy.report` module with custom HTML plots, native [AG-Grid](https://www.ag-grid.com/) implementation, stat boxes, etc. 
!!!

## Using Harpy reports
Harpy workflows still create reports during workflows, as they did before, but they aren't rendered as HTML documents
like they were before. You could open the `.ipynb` files in Jupyter/VScode/etc., but to really get the benefits of the
reports, you need to lean on `jupyter-book` to compile/render everything into a single website. That can be done with
`harpy report`, which, when run in the project root directory, will create the necessary MyST configurations and download the website
templates/assets if necessary, then start a local liveserver to view the reports.

```bash
harpy report <options> DIRECTORY
```

{.compact}
| argument    |   default   | description |
|:-----------------------|:--------------:|:--------------------------------------------------------|
| `DIRECTORY` | `.` | Path of where to setup configs and launch server |
| `--debug` `-d` | | Dump all of jupyterbook's output to the terminal |
| `--headless` `-h` | |  Run the server in headless mode, with only the content server started |
| `--md` `-m` | | Also scan for markdown files (`.md`) |
| `--clear-cache` `-c` |  False | Remove `_build` directory prior to server launch |
| `--port` `-p` | | Run the application server from the specified port number |
| `--refresh` `-r` | 0 | Refresh interval, in seconds (disabled with `0`) |
| `--server-port` `-s` | | Run the content server from the specified port number |

## Automate a report website
The report setup and system lends itself well to easily build [a persistent report website](https://pdimens.github.io/RS-degradation/) via GitHub Pages.
This requires one additional file: a GitHub Actions workflow that is triggered on push events (or whatever
you configure it to). This can be added using:
```bash
harpy template report --action
```
The `template report` command will make sure all the necessary files are there, which include:
- `myst.yml`
- `.report/index.md`
- `.report/favicon.ico`
- `.report/logo.png`
- `.github/workflows/buildreports.yml`
- `.github/dependabot.yml` 
