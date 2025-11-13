# Snakemake Logger Plugin: Rich

**Note:** This plugin is currently in development and may not be fully stable. Use at your own discretion and report any issues to the repository's issue tracker.

## Introduction

A logging plugin for Snakemake that utilizes [rich](https://github.com/Textualize/rich) for enhanced terminal styling and progress bars. 

## Usage

1. Install via pip:
   ```
   pip install snakemake-logger-plugin-rich
   ```

2. Run Snakemake with the `--logger rich` option:
   ```
   snakemake --logger rich
   ```

## Development

This project uses [pixi](https://github.com/prefix-dev/pixi) for environment management.

### Setting up the development environment

1. Fork the repository

2. Clone the repository:
   ```
   git clone https://github.com/<your-username>/snakemake-logger-plugin-rich.git
   cd snakemake-logger-plugin-rich
   ```

3. Install dependencies using pixi:
   ```
   pixi install
   ```

4. Activate the development environment:
   ```
   pixi shell -e dev
   ```

### Available development tasks

Run these commands with `pixi run`:

- **Demo**: `pixi run demo [snakefile]` - Runs a complete demo workflow and cleans up afterward. Uses `demo/Snakefile` by default, or specify a custom Snakefile path
- **Run Demo**: `pixi run run-demo [snakefile]` - Runs a demonstration Snakemake workflow using the plugin. Uses `demo/Snakefile` by default, or specify a custom Snakefile path
- **Dry Run**: `pixi run dryrun [snakefile]` - Performs a dry run of the demo workflow. Uses `demo/Snakefile` by default, or specify a custom Snakefile path
- **Clean Demo**: `pixi run clean-demo` - Cleans up demo output files
- **Quality Control**: `pixi run qc` - Runs formatting, linting, and type checking
- **Format**: `pixi run format` - Format code with ruff
- **Lint**: `pixi run lint` - Lint code with ruff
- **Lint Fix**: `pixi run lint-fix` - Lint and auto-fix code with ruff
- **Type Check**: `pixi run type-check` - Type check with mypy

### Testing the plugin

To test the plugin with the demo workflow:

```
pixi run demo
```

Or to run just the demo workflow without cleanup:

```
pixi run run-demo
```

To test with a custom Snakefile:

```
pixi run demo path/to/your/Snakefile
```
