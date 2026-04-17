{{ PYTHON }} -m pip install . --no-deps --no-build-isolation --no-cache-dir -vvv

## build Go binaries
{
    cd harpy/utils
    go build -C stagger -o ../gih-stagger -ldflags='-s -w' stagger.go
    go build -C convert -o ../gih-convert -ldflags='-s -w' convert.go
    go build -C standardize -o ../djinn-standardize -ldflags='-s -w' standardize.go 
    chmod +x gih-stagger gih-convert djinn-standardize
    mv gih-stagger gih-convert djinn-standardize ${PREFIX}/bin/
}

## CLI completions
mkdir -p "$PREFIX/etc/conda/activate.d"
mkdir -p "$PREFIX/share/harpy"

# Pre-generate completion scripts
_HARPY_COMPLETE=bash_source harpy > "$PREFIX/share/harpy/complete.bash"
_HARPY_COMPLETE=zsh_source harpy  > "$PREFIX/share/harpy/complete.zsh"
_HARPY_COMPLETE=fish_source harpy > "$PREFIX/share/harpy/complete.fish"

# Write the activate hook
cat > "$PREFIX/etc/conda/activate.d/harpy-completion.sh" << 'EOF'
_shell=$(ps -p $$ -o comm= 2>/dev/null || basename "$SHELL")
case "$_shell" in
  bash) source "$CONDA_PREFIX/share/harpy/complete.bash" ;;
  zsh)  source "$CONDA_PREFIX/share/harpy/complete.zsh"  ;;
  fish) source "$CONDA_PREFIX/share/harpy/complete.fish" ;;
esac
EOF