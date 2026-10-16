# Sourced on activation of the harpy environment (conda activate.d).
# Enables tab-completion for harpy in fish from the script installed under $CONDA_PREFIX/share.
set -l _harpy_completions "$CONDA_PREFIX/share/fish/vendor_completions.d"
if test -d "$_harpy_completions"; and not contains -- "$_harpy_completions" $fish_complete_path
    set -a fish_complete_path "$_harpy_completions"
end
