{{ PYTHON }} -m pip install . --no-deps --no-build-isolation --no-cache-dir -vvv

python -m ipykernel install --user

## preproc commands
{
    cd harpy/utils/stagger
    go mod tidy && go build -ldflags="-s -w" -o ../gih-stagger stagger.go
    cd ../convert
    go mod tidy && go build -ldflags="-s -w" -o ../gih-convert convert.go
    cd .. && chmod +x gih-stagger gih-convert
    mv gih-stagger gih-convert ${CONDA_PREFIX}/bin/
}
