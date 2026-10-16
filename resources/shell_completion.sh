# Sourced on activation of the harpy environment (conda activate.d / pixi activation script).
# Enables tab-completion for harpy, harpy-utils, and hv (bash, zsh, fish) from the scripts installed under $CONDA_PREFIX/share.
#
# Some tools (`pixi shell`) run activation scripts in a subprocess and only keep the environment
# variables they set, so code that registers completions here can't be relied upon. Instead,
# point the shells' own completion lookup at this environment via an environment variable:
#   XDG_DATA_DIRS: bash-completion (bash) and fish look for completions in <dir>/... for each entry
# Do NOT set FPATH for zsh: it isn't needed (pixi adds share/zsh/site-functions to fpath and runs compinit
# by itself), and an exported FPATH replaces zsh's whole function path (including what ~/.zshrc and
# frameworks like oh-my-zsh add), which breaks the shell.
_harpy_share="${CONDA_PREFIX-}/share"
if [ -d "${_harpy_share}" ]; then
  case ":${XDG_DATA_DIRS-}:" in
    *":${_harpy_share}:"*) ;;
    *) export XDG_DATA_DIRS="${XDG_DATA_DIRS:-/usr/local/share:/usr/share}:${_harpy_share}" ;;
  esac
fi

# When this does run inside the user's own interactive shell (conda), also register completion
# right away, which additionally works without bash-completion in bash.
case $- in
  *i*)
    if [ -n "${BASH_VERSION-}" ]; then
      for _harpy_program in harpy harpy-utils hv; do
        if [ -f "${_harpy_share}/bash-completion/completions/${_harpy_program}" ]; then
          . "${_harpy_share}/bash-completion/completions/${_harpy_program}"
        fi
      done
      unset _harpy_program
    elif [ -n "${ZSH_VERSION-}" ]; then
      # compinit (and therefore compdef) normally isn't set up until ~/.zshrc runs, after this,
      # so register right before the first prompt instead.
      _harpy_completion_zsh() {
        add-zsh-hook -d precmd _harpy_completion_zsh
        unset -f _harpy_completion_zsh
        if type compdef >/dev/null 2>&1; then
          for _harpy_program in harpy harpy-utils hv; do
            if [ -f "${CONDA_PREFIX-}/share/zsh/site-functions/_${_harpy_program}" ]; then
              . "${CONDA_PREFIX}/share/zsh/site-functions/_${_harpy_program}"
            fi
          done
          unset _harpy_program
        fi
      }
      autoload -Uz add-zsh-hook && add-zsh-hook precmd _harpy_completion_zsh
    fi
    ;;
esac
unset _harpy_share
