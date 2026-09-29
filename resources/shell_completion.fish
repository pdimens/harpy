# Sourced on activation of the harpy environment (conda activate.d).
# Enables tab-completion for harpy in fish using the script that is written to
# $CONDA_PREFIX/share/harpy/ when harpy is built (see generate_completion.sh).
if status is-interactive; and test -f "$CONDA_PREFIX/share/harpy/complete.fish"
    source "$CONDA_PREFIX/share/harpy/complete.fish"
end
