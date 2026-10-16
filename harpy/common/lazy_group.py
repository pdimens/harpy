"""Click groups that import their subcommands only when they are needed.

Every Harpy command module pulls in a fair amount of machinery (pysam, yaml, rich, ...),
so importing all of them just to run `harpy --version`, print the top-level help, or
tab-complete a subcommand name is wasted work. `LazyGroup` stores an import path per
subcommand and only imports the module of the subcommand that is actually being run.

This module must stay light: it is imported by the entry points and should only depend on click.
For rich-click help formatting, combine it with `rich_click.RichGroup` at the call site.
"""

import importlib
from typing import NamedTuple

import click


class LazySpec(NamedTuple):
    """Where to find a subcommand and how to describe it without importing it."""

    target: str
    """Import path of the command in the form `package.module:attribute`"""
    help: str | None = None
    """
    First paragraph of the command's help text. When given, the group can list the
    command (`--help`, shell completion) without importing it. When `None`, listing the
    command imports it.
    """
    hidden: bool = False
    """Mirror of the command's `hidden` flag, only used when `help` is given."""


class LazyCommand(click.Command):
    """
    A lightweight stand-in for a command that hasn't been imported yet. It carries just enough
    information (name, help, hidden) to be listed in a group's help or shell completion.
    It is never invoked: `LazyGroup.resolve_command` swaps it for the real command.
    """

    def __init__(self, name: str, spec: LazySpec):
        super().__init__(name, help=spec.help, hidden=spec.hidden)
        self.spec = spec

    def load(self) -> click.Command:
        """Import and return the real command."""
        module_path, attribute = self.spec.target.split(":")
        return getattr(importlib.import_module(module_path), attribute)


class LazyGroup(click.Group):
    """
    A click group whose subcommands are declared as `{name: LazySpec}` and imported on demand.

    - `list_commands` never imports anything
    - `get_command` returns a `LazyCommand` stand-in for specs that define `help`, otherwise the real command
    - `resolve_command` (what click uses to dispatch and to descend during shell completion) always
      returns the real command

    Anything registered the usual way with `add_command` keeps working alongside the lazy ones.
    """

    def __init__(self, *args, lazy_commands: dict[str, LazySpec] | None = None, **kwargs):
        super().__init__(*args, **kwargs)
        self.lazy_commands: dict[str, LazySpec] = lazy_commands or {}
        self._lazy_cache: dict[str, click.Command] = {}

    def list_commands(self, ctx: click.Context) -> list[str]:
        return sorted({*super().list_commands(ctx), *self.lazy_commands})

    def get_command(self, ctx: click.Context, cmd_name: str) -> click.Command | None:
        if cmd_name not in self.lazy_commands:
            return super().get_command(ctx, cmd_name)
        if cmd_name not in self._lazy_cache:
            stub = LazyCommand(cmd_name, self.lazy_commands[cmd_name])
            self._lazy_cache[cmd_name] = stub if stub.spec.help is not None else stub.load()
        return self._lazy_cache[cmd_name]

    def resolve_command(self, ctx: click.Context, args: list[str]):
        cmd_name, cmd, remaining = super().resolve_command(ctx, args)
        if isinstance(cmd, LazyCommand):
            cmd = self._lazy_cache[cmd.name] = cmd.load()
        return cmd_name, cmd, remaining
