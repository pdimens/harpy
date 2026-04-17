_shell=$(ps -p $$ -o comm= 2>/dev/null || basename "$SHELL")
case "$_shell" in
  bash) source "$CONDA_PREFIX/share/harpy/complete.bash" ;;
  zsh)  source "$CONDA_PREFIX/share/harpy/complete.zsh"  ;;
  fish) source "$CONDA_PREFIX/share/harpy/complete.fish" ;;
esac