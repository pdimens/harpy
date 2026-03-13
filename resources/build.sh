{{ PYTHON }} -m pip install . --no-deps --no-build-isolation --no-cache-dir -vvv

python -m ipykernel install --user

## preproc commands
{
    cd harpy/utils
    go mod tidy
    go build -ldflags="-s -w" -o gih-stagger stagger/stagger.go
    go build -ldflags="-s -w" -o gih-convert convert/convert.go
    chmod +x gih-stagger gih-convert
    cp gih-stagger gih-convert ${PREFIX}/bin/
}
