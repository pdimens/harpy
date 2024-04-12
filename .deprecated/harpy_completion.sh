#! /usr/bin/env bash

shell=$(basename $SHELL)

if [ "$shell" = "bash" ] || [ "$shell" = "sh" ]; then
    eval "$(_HARPY_COMPLETE=bash_source harpy)"
elif [ "$shell" = "zsh" ]; then
    eval "$(_HARPY_COMPLETE=zsh_source harpy)"
elif [ "$shell" = "fish" ]; then
    _HARPY_COMPLETE=fish_source harpy | source
else
    echo -n ""
fi