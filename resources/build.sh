{{ PYTHON }} -m pip install . --no-deps --no-build-isolation --no-cache-dir -vvv

## shell completion scripts (bash, zsh, fish) for harpy, harpy-utils, and hv, installed where each shell looks for them
mkdir -p ${PREFIX}/share/bash-completion/completions ${PREFIX}/share/zsh/site-functions ${PREFIX}/share/fish/vendor_completions.d
for _program in harpy harpy-utils hv; do
    {{ PYTHON }} -m harpy completion bash ${_program} > ${PREFIX}/share/bash-completion/completions/${_program}
    {{ PYTHON }} -m harpy completion zsh ${_program} > ${PREFIX}/share/zsh/site-functions/_${_program}
    {{ PYTHON }} -m harpy completion fish ${_program} > ${PREFIX}/share/fish/vendor_completions.d/${_program}.fish
done

## build Go binaries
{
    cd harpy/utils
    go build -C stagger -o ../gih-stagger -ldflags='-s -w' stagger.go
    go build -C convert -o ../gih-convert -ldflags='-s -w' convert.go
    go build -C standardize -o ../djinn-standardize -ldflags='-s -w' standardize.go 
    chmod +x gih-stagger gih-convert djinn-standardize
    mv gih-stagger gih-convert djinn-standardize ${PREFIX}/bin/
}

## activate/deactive processes
mkdir -p $PREFIX/etc/conda/activate.d
mkdir -p $PREFIX/etc/conda/deactivate.d

## Keep these two hooks identical to resources/shell_completion.{sh,fish} (checked by test/cli/check_lazy_cli.py)
cat > ${PREFIX}/etc/conda/activate.d/harpy-completion.sh <<'EOF_COMPLETION_SH'
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
EOF_COMPLETION_SH

cat > ${PREFIX}/etc/conda/activate.d/harpy-completion.fish <<'EOF_COMPLETION_FISH'
# Sourced on activation of the harpy environment (conda activate.d).
# Enables tab-completion for harpy in fish from the script installed under $CONDA_PREFIX/share.
set -l _harpy_completions "$CONDA_PREFIX/share/fish/vendor_completions.d"
if test -d "$_harpy_completions"; and not contains -- "$_harpy_completions" $fish_complete_path
    set -a fish_complete_path "$_harpy_completions"
end
EOF_COMPLETION_FISH

cat > ${PREFIX}/etc/conda/activate.d/harpy-activate.sh <<'EOF'  
export _HARPY_OLD_JUPYTER_NOTARY_DB="${JUPYTER_NOTARY_DB-__UNSET__}"  
export JUPYTER_NOTARY_DB=':memory:'

python -m ipykernel install --prefix "$CONDA_PREFIX" --name ipython-harpy \
    --display-name "Python (harpy)"
EOF
  
cat > ${PREFIX}/etc/conda/deactivate.d/harpy-deactivate.sh <<'EOF'  
if [ "${_HARPY_OLD_JUPYTER_NOTARY_DB-__UNSET__}" = "__UNSET__" ]; then  
  unset JUPYTER_NOTARY_DB  
else  
  export JUPYTER_NOTARY_DB="${_HARPY_OLD_JUPYTER_NOTARY_DB}"  
fi  
unset _HARPY_OLD_JUPYTER_NOTARY_DB  
EOF
