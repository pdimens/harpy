import base64
import gzip
import json
from IPython.display import display, HTML

class JSFunction:
    """Wraps a raw JS string so it can be injected without JSON quoting."""
    def __init__(self, js):
        """
        Store a JavaScript snippet for raw insertion into serialized column definitions.
        
        Parameters:
        	js (str): JavaScript code to embed; leading and trailing whitespace will be removed and the result stored on `self.js`.
        """
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
        """
        Initialize ITable configuration and prepare column definitions and row data for rendering.
        
        Parameters:
            df: A DataFrame-like object (pandas or polars) providing column names and row data.
            filename (str): Default filename used by the grid export action.
            fixedcols (int): Number of leftmost columns to pin in the grid.
            compress (bool): If True and row data is non-empty, prepare a gzip+base64 compressed columnar payload for embedding; otherwise embed raw row records.
        
        Notes:
            - Builds `col_defs` from `df.columns`, pinning the first `fixedcols` columns.
            - Normalizes input data to a list of row dictionaries.
            - Assigns unique `grid_id` and `grid_ref` identifiers derived from the DataFrame object id.
            - When compression is enabled and data exists, stores the compressed payload in `_compressed_payload` and leaves `row_data` as None; otherwise stores the raw row list in `row_data`.
        """
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
        Convert a list of row-oriented records into a gzip-compressed, base64-encoded columnar JSON payload.
        
        The produced JSON object has the shape {"cols": [col1, col2, ...], "data": {col1: [v1, v2, ...], ...}}; it is UTF-8 encoded, gzip-compressed (level 9), and base64-encoded for safe embedding in HTML/JS.
        
        Returns:
            A base64-encoded ASCII string containing the gzip-compressed JSON payload.
        """
        cols = list(raw[0].keys())
        columnar = {col: [row[col] for row in raw] for col in cols}
        payload = json.dumps(
            {"cols": cols, "data": columnar}, default=str
        ).encode("utf-8")
        compressed = gzip.compress(payload, compresslevel=9)
        return base64.b64encode(compressed).decode("ascii")


    def _serialize_col_defs(self) -> str:
        """
        Serialize the table's column definition objects into a JavaScript array literal, embedding JSFunction values as raw JavaScript.
        
        Returns:
            str: A JavaScript array literal string representing the column definitions; values of type `JSFunction` are inserted without JSON quoting while all other values are JSON-encoded.
        """
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
        Produce a JavaScript snippet that defines the grid's `rowData`.
        
        When a compressed payload is present the snippet performs base64 decoding and gzip decompression in the browser, reconstructs row-oriented records from a columnar structure, and assigns the resulting array to `rowData`. When no compressed payload is used the snippet is a `const rowData = ...;` assignment containing the raw JSON-serialized rows.
        
        Returns:
            str: JavaScript code that, when executed, defines the `rowData` variable for AG Grid.
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
        """
        Generate the AG Grid HTML for the table and either return it or display it in the current IPython output.
        
        Parameters:
            html (bool): If True, return the generated HTML string; if False, display the HTML in the current IPython output.
        
        Returns:
            The generated HTML string when `html` is True, `None` otherwise.
        """

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

