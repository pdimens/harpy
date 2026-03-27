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

