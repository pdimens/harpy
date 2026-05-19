import base64
import gzip
import json
from IPython.display import display, HTML

class JSFunction:
    """Wraps a raw JS string so it can be injected without JSON quoting."""
    def __init__(self, js):
        self.js = js.strip()

class ITable:
    """
    Render a pandas/polars DataFrame as an AG Grid table using self-contained HTML.
    Works in Jupyter notebooks and static MyST/Jupyter Book pages.

    Parameters
    ----------
    df          : pandas or polars DataFrame to display
    filename    : default filename used when exporting to CSV
    fixedcols   : number of columns to pin/freeze on the left
    compress    : if True (default), embed row data as gzip+base64 columnar JSON
                  instead of raw JSON — typically smaller file size
    """

    def __init__(self, df, filename: str, fixedcols: int = 0, compress: bool = True):
        self.theme: str = "ag-theme-quartz"
        self.row_height: int = 28
        self.header_height: int = 32
        self.filename: str = filename
        self.compress: bool = compress

        self.col_defs = [
            {"field": col, **({"pinned": "left"} if i < fixedcols else {})}
            for i, col in enumerate(df.columns)
        ]

        # Normalise to list-of-dicts regardless of pandas/polars
        raw = (
            df.to_dicts()
            if hasattr(df, "to_dicts")
            else df.to_dict(orient="records")
        )

        self.grid_id = f"grid-{id(df)}"
        self.grid_ref = f"aggrid_{id(df)}"

        if compress and raw:
            self._compressed_payload = self._build_compressed_payload(raw)
            self.row_data = None          # not used in compressed path
        else:
            self._compressed_payload = None
            self.row_data = raw

    @staticmethod
    def _build_compressed_payload(raw: list[dict]) -> str:
        """
        1. Repack row-oriented data into a columnar dict  →  far fewer repeated keys
        2. JSON-serialise the columnar structure
        3. gzip compress (level 9)
        4. base64-encode so it's safe to embed in a <script> string

        Returns the base64 string.
        """
        cols = list(raw[0].keys())
        columnar = {col: [row[col] for row in raw] for col in cols}
        payload = json.dumps(
            {"cols": cols, "data": columnar}, default=str
        ).encode("utf-8")
        compressed = gzip.compress(payload, compresslevel=9)
        return base64.b64encode(compressed).decode("ascii")


    def _serialize_col_defs(self) -> str:
        '''
        Column-def serialisation (preserves raw JSFunction values)
        '''
        parts = []
        for col in self.col_defs:
            field_parts = []
            for k, v in col.items():
                if isinstance(v, JSFunction):
                    field_parts.append(f'"{k}": {v.js}')
                else:
                    field_parts.append(f'"{k}": {json.dumps(v)}')
            parts.append("{" + ", ".join(field_parts) + "}")
        return "[" + ", ".join(parts) + "]"

    def _row_data_js(self) -> str:
        """
        Returns a JS snippet that declares `rowData` as an array of objects.

        Compressed path  → async decompress from embedded base64/gzip string.
        Uncompressed path → plain JSON literal
        """
        if self._compressed_payload:
            return f"""
                // ---------- decompress columnar gzip payload ----------
                const _b64 = "{self._compressed_payload}";
                const _binary = Uint8Array.from(atob(_b64), c => c.charCodeAt(0));
                const _ds = new DecompressionStream("gzip");
                const _writer = _ds.writable.getWriter();
                _writer.write(_binary);
                _writer.close();
                const _buf = await new Response(_ds.readable).arrayBuffer();
                const _pkg = JSON.parse(new TextDecoder().decode(_buf));

                // Reconstruct row-oriented array from columnar structure
                const _cols = _pkg.cols;
                const _data = _pkg.data;
                const rowData = _data[_cols[0]].map((_, i) =>
                    Object.fromEntries(_cols.map(c => [c, _data[c][i]]))
                );
                // ------------------------------------------------------
            """
        # Fallback: raw JSON
        return f"const rowData = {json.dumps(self.row_data, default=str)};"

    def render(self, html: bool = False):
        """Create the AG-Grid HTML and render it (or return the HTML string)."""

        _html = f"""
        <style>
            .ag-theme-quartz, .ag-theme-quartz-dark {{
                --ag-font-family: sans-serif;
            }}
        </style>

        <button
            onclick="(function(){{ var g = window.{self.grid_ref}; if(g) g.exportDataAsCsv({{suppressQuotes: true, fileName: '{self.filename}'}}); }})()"
            style="margin-bottom: 8px; padding: 4px 12px; cursor: pointer;"
        >
            Export CSV
        </button>

        <div id="{self.grid_id}" class="{self.theme}" style="width: 100%; overflow-x: auto"></div>

        <script>
        (async function () {{

                    /* ── dynamic loader: awaitable, idempotent ── */
            await new Promise((resolve, reject) => {{
                if (typeof agGrid !== "undefined") {{ resolve(); return; }}
                const s = document.createElement("script");
                s.src = "https://cdn.jsdelivr.net/npm/ag-grid-community/dist/ag-grid-community.min.js";
                s.onload = resolve;
                s.onerror = reject;
                document.head.appendChild(s);
            }});

            {self._row_data_js()}

            const columnDefs = {self._serialize_col_defs()};

            const gridOptions = {{
                columnDefs: columnDefs,
                rowData: rowData,
                defaultColDef: {{
                    sortable: true,
                    filter: true,
                    resizable: true,
                }},
                autoSizeStrategy: {{
                    type: "fitCellContents",
                    defaultMaxWidth: 170,
                    defaultMinWidth: 90,
                }},
                domLayout: "autoHeight",
                animateRows: false,
                pagination: true,
                paginationPageSize: 20,
                rowHeight: {self.row_height},
                headerHeight: {self.header_height},
            }};

            const container = document.getElementById("{self.grid_id}");

            function syncTheme() {{
                container.setAttribute(
                    "data-ag-theme-mode",
                    document.documentElement.classList.contains("dark") ? "dark-blue" : "light"
                );
            }}

            requestAnimationFrame(() => {{
                requestAnimationFrame(() => {{
                    syncTheme();
                    window.{self.grid_ref} = agGrid.createGrid(container, gridOptions);
                    new MutationObserver(syncTheme).observe(
                        document.documentElement,
                        {{ attributes: true, attributeFilter: ["class"] }}
                    );
                }});
            }});

        }})();
        </script>
        """

        if html:
            return _html
        return display(HTML(_html))

