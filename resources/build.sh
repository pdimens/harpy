{{ PYTHON }} -m pip install . --no-deps --no-build-isolation --no-cache-dir -vvv

## stagger command
cd stagger && go build -ldflags="-s -w" -o stagger-gih stagger.go
chmod +x stagger-gih && cp stagger-gih TO SOMEWHWERE