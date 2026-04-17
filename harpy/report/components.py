import base64
import html
from pathlib import Path
from datetime import datetime
from io import BytesIO
from IPython.display import display, HTML, Image
from PIL import Image as PImage
import uuid
import altair as alt
import pandas as pd
from .theme import palette
import json

class StatsBox:
    '''
    The boxes class that draws boxes at the top of the reports.
    '''
    def __init__(self, height: int = 85, label_fontsize=14, value_fontsize=36, text_gap: int = 1):
        self.boxes: list[str] = []
        self.height = height
        self.label_fontsize = label_fontsize
        self.value_fontsize = value_fontsize
        self.pm_fontsize = value_fontsize // 3
        self.text_gap = text_gap
        self.html = (
            '<div style="background-color: transparent ; '
            f'height: {self.height}px; display: inline-flex; flex-direction: column; white-space: nowrap; '
            'align-items: flex-end; justify-content: center; border-radius: 5px; color: #000000; '
            'padding: 15px; box-sizing: border-box; text-align: right;">'
            f'<div style="font-size: {self.label_fontsize}px; font-weight: 400;' 
            f'margin-bottom: {self.text_gap}px; color: var(--tw-prose-body, #666666);">{{label}}</div>'
            f'<div style="font-size: {self.value_fontsize}px; font-weight: 400; color: {{color}};">{{value}}'
            f'<span style="font-size: {self.value_fontsize // 2}px; color: {{color}};">{{plus_minus}}</span>'
            '</div>'
            '</div>'
        )

    def add(self, value, label, plus_minus: str|int|float = '', textcol=None):
        '''
        Return the html of a colored box object with `value` and `label`. 
        
        Args:
            value: Main metric to display
            label: Title of the metric
            plus_minus: smaller text to add if there is a plus-minus value
            textcol: color to make the value text
        '''
        _val = value if isinstance(value, str) else f"{value:,.2f}".rstrip('0').rstrip('.')
        _col = textcol or 'var(--tw-prose-body, #666666)'
        pm_str = f'±{plus_minus}' if plus_minus else ''
        self.boxes.append(
            self.html.format(label = label, value = _val, color = _col, plus_minus=pm_str)
        )
        return self

    def conditional(self, value, label, cutoff: int|float, lower_bad: bool = True, as_percent:bool = False, plus_minus: str|int|float = ''):
      '''
      Return the html of a colored box object with `value` and `label`. Use `as_percent` to multiply
      the value by 100 for printing purposes.
      The `color` is either yellow or green depending on what is determined better or worse than the `cutoff`:
      - `lower_bad=True`: `color` = yellow when value < cutoff (default)
      - `lower_bad=False`: `color` = yellow when value >= cutoff 
      '''
      if lower_bad:
        color = "#f6ab3c" if value < cutoff else None
      else:
        color = "#f6ab3c" if value >= cutoff else None
      return self.add(value if not as_percent else f"{value * 100}%", label, plus_minus=plus_minus, textcol = color)

    def render(self, gap: int = 5):
        '''Display all boxes in a horizontal row
        
        Args:
            gap: Space between boxes in pixels (default: 5)
        '''
        container_html = f'''<div style="display: flex; gap: {gap}px; 
            flex-wrap: wrap;">{" ".join(self.boxes)}</div>'''
        display(HTML(container_html))

def print_time(*args):
    '''HTML-print the time and all arguments with line breaks between them'''
    _now = datetime.now().strftime("🗓️ %d %B, %Y 🕔 %H:%M")
    _html = "<p>{}<p>".format("<br>".join([_now] + [str(i) for i in [*args]]))
    return display(HTML(_html))

def print_html(*args):
    '''HTML-print all arguments with line breaks between them'''
    display(HTML(
        "<p>{}</p>".format("<br>".join(str(i) for i in [*args]))
    ))

def embed_image(x: str, scale: float = 0.5):
    '''Rescale an PNG image and embed it into the notebook'''
    image = PImage.open(x)
    if scale < 1:  
        image = image.resize((int(image.size[0] * scale), int(image.size[1] * scale)), PImage.LANCZOS)  
    buf = BytesIO()  
    image.save(buf, format = "PNG" if "png" in x.lower() else "jpeg")  
    return Image(data=buf.getvalue(), format="png" if "png" in x.lower() else "jpeg", embed=True)

def image_viewer(label: str, image_dir: str, pattern: str, sortkey = None, thing_to_select: str = "sample", recursive: bool = False, scale: float = 1.0, option_key = None):  
    '''
    Create a javascript image viewer with a file picker.
    Images are embedded (thus can safely have their original files deleted) and can be down-scaled.
    '''
    imgfmt = "image/png" if "png" in pattern else "image/jpg"
    options_parts = []
    paths = Path(image_dir).glob(pattern) if not recursive else Path(image_dir).rglob(pattern)  
    paths = sorted(paths) if not sortkey else sorted(paths, key = sortkey)
    uid = uuid.uuid4().hex[:8]
    images_dict = {}
    option_key = option_key or (lambda p: p.stem if not recursive else p.parents[1].name) 
    for p in paths:
        pname = option_key(p)
        safe_value = html.escape(pname, quote=True)  
        safe_label = html.escape(pname)  
        options_parts.append(f'<option value="{safe_value}">{safe_label}</option>') 
        #options_parts.append(f'<option value="{pname}">{pname}</option>')
        image = PImage.open(p)
        if scale < 1:
            image = image.resize((int(image.size[0] * scale), int(image.size[1] * scale)), PImage.LANCZOS)
        buf = BytesIO()
        image.save(buf, format="PNG" if "png" in pattern else "jpeg")
        buf.seek(0)
        encoded = base64.b64encode(buf.read()).decode()
        images_dict[pname] = f"data:{imgfmt};base64,{encoded}"
    images_js = json.dumps(images_dict)
    options_html = "\n".join(options_parts)

    _html = f"""<!DOCTYPE html>
        <html>
        <head>
        <style>
        body {{ font-family: sans-serif; padding: 1rem; }}
        #viewer-{uid} {{ margin-bottom: 1rem; }}
        #viewer-{uid} img {{ max-width: 100%; }}
        #selector-{uid} {{ display: flex; align-items: center; gap: 0.75rem; margin-bottom: 1rem; }}
        #selector-{uid} label {{ font-size: 0.9rem; font-weight: bold; }}
        #file-select-{uid} {{
            padding: 6px 10px;
            font-size: 0.9rem;
            border: 1px solid #ccc;
            border-radius: 4px;
            cursor: pointer;
            min-width: 200px;
        }}
        </style>
        </head>
        <body>
        <div id="selector-{uid}">
        <label for="file-select-{uid}">{html.escape(label)} </label>  
        <select id="file-select-{uid}" onchange="showImage_{uid}(this.value)">
            <option value="" disabled selected>Select a {html.escape(thing_to_select)}</option> 
            {options_html}
        </select>
        </div>
        <div id="viewer-{uid}"><img id="display-{uid}" src="" alt="Select a {thing_to_select}" /></div>
        <script>
        const images_{uid} = {images_js};

        function showImage_{uid}(name) {{
            document.getElementById("display-{uid}").src = images_{uid}[name];
        }}
        </script>
        </body>
        </html>"""
    display(HTML(_html))

#def standard_itable(data, filename:str, caption: str|None= None, fixedcols: int|None = None, coldefs = [], html: bool = False):
#    '''
#    Boilerplate version of `itables.show` to reduce in-report boilerplate. Use `filename` to specify the output prefix for
#    exports, `caption` to set an optional caption, `fixcols` to freeze this-many columns on the left, and `coldefs` to
#    populate `columnDefs` with wacky JS.
#    '''
#    return itables.show(
#        data,
#        caption = caption,
#        buttons = [
#            "pageLength", "colvis",
#            {"extend": "csvHtml5", "title": filename},
#            {"extend": "excelHtml5", "title": filename}
#        ],
#        fixedColumns = {"left": fixedcols} if fixedcols else {},
#        ordering = {"indicators": False, "handler": False},
#        columnControl = [["order", "searchDropdown"]],
#        showIndex= False,
#        scrollX =  True,
#        autoWidth=True,
#        searching = False,
#        style = "width:100%",
#        classes = "display nowrap compact",
#        allow_html=html,
#        columnDefs = coldefs
#    )
#

class JSFunction:
    """Wraps a raw JS string so it can be injected without JSON quoting."""
    def __init__(self, js):
        self.js = js.strip()

class ITable:
    """
    Render a pandas/polars DataFrame as an AG Grid table using self-contained HTML.
    Works in Jupyter notebooks and static MyST/Jupyter Book pages. `df` is the DataFrame
    to display, and `frozen_cols` is the number of columns to pin/freeze on the left. Includes
    additional configurable fields `theme`, `row_height`, and `header_height`.
    """
    def __init__(self, df, filename: str, fixedcols: int = 0):
        self.theme: str = "ag-theme-quartz"
        self.row_height: int = 28
        self.header_height: int = 32
        self.filename: str = filename
        self.col_defs = [{"field": col, **({"pinned": "left"} if i < fixedcols else {})} for i, col in enumerate(df.columns)]
        self.row_data = df.to_dicts() if hasattr(df, "to_dicts") else df.to_dict(orient="records")
        self.grid_id = f"grid-{id(df)}"
        self.grid_ref = f"aggrid_{id(df)}"
        # save autosize coldef for a rainy day?
        ##autosize = {"autoSizeStrategy": "AutoSizeStrategy" = {"type": "fitCellContents", "defaultMaxWidth": 150, "defaultMinWidth": 80, "scaleUpToFitGridWidth": true}}
        
    def _serialize_col_defs(self):
        """Serialize col_defs, injecting JSFunction values as raw JS."""
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

    def render(self, html: bool = False):
        '''Create the AG-Grid HTML and render it'''
        _html = f"""
        <button
            onclick="(function(){{ var g = document.{self.grid_ref}; if(g) g.exportDataAsCsv({{suppressQuotes: true, fileName : "{self.filename}"}}); }})()"
            style="margin-bottom: 8px; padding: 4px 12px; cursor: pointer;"
        >
            Export CSV
        </button>

        <div id="{self.grid_id}" class="{self.theme}" style="width: 100%; overflow-x: auto"></div>
        <script src="https://cdn.jsdelivr.net/npm/ag-grid-community/dist/ag-grid-community.min.js"></script>

        <script>
            (function() {{
                const columnDefs = {self._serialize_col_defs()};
                const rowData = {json.dumps(self.row_data, default=str)};

                const gridOptions = {{
                    columnDefs: columnDefs,
                    rowData: rowData,
                    defaultColDef: {{
                        sortable: true,
                        filter: true,
                        resizable: true,
                    }},
                    autoSizeStrategy: AutoSizeStrategy = {{
                        type: "fitCellContents",
                        defaultMaxWidth: 150,
                        defaultMinWidth: 80,
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
        
                function numberParser(params) {{
                const newValue = params.newValue;
                let valueAsNumber;
                if (newValue === null || newValue === undefined || newValue === "") {{
                    valueAsNumber = null;
                }} else {{
                    valueAsNumber = parseFloat(params.newValue);
                }}
                return valueAsNumber;
                }}

                function tryInit(attempts) {{
                    if (typeof agGrid !== "undefined") {{
                        requestAnimationFrame(() => {{
                            requestAnimationFrame(() => {{
                                syncTheme();
                                document.{self.grid_ref} = agGrid.createGrid(container, gridOptions);
                                new MutationObserver(syncTheme).observe(document.documentElement, {{
                                    attributes: true, attributeFilter: ["class"]
                                }});
                            }});
                        }});
                    }} else if (attempts > 0) {{
                        setTimeout(() => tryInit(attempts - 1), 100);
                    }}
                }}
                tryInit(20);
            }})();
        </script>    """
        if html:
            return _html
        return display(HTML(_html))


def piechart(df: pd.DataFrame, title: str = "", goodval: str|None = None, lgtitle: bool = False, lgpos: str = "bottom"):
    total = df['count'].sum()
    _names = list(df.columns)
    _n = set(df[_names[0]])
    _color = alt.Color(
        f'{_names[0]}:N', title = _names[0].title(),
        legend=alt.Legend(orient=lgpos, title = None) if not lgtitle else alt.Legend(orient=lgpos)
        )
    if goodval:
        badval = [i for i in _n if i != goodval][0]
        _color = _color.scale(domain=[goodval, badval], range= ["#68ae6b","#C9C9C9"])
    return (
    alt.Chart(df)
    .transform_calculate('percent', f'datum.count / {total}')
    .mark_arc()
    .encode(
        theta=alt.Theta('count:Q'),
        color = _color,
        tooltip = [
            alt.Tooltip(f'{_names[0]}:O', title = _names[0].title()),
            alt.Tooltip('percent:Q', title = "Percent", format = ".2%"),
            alt.Tooltip("count:Q", title = "Count")
        ]
    )
    .properties(
        height = 100, width = 150,
        title = alt.Title(title, align = 'center', anchor = 'middle')
    )
)

def convolutionpie(df: pd.DataFrame, title: str = "", lgpos: str = 'bottom'):
    total = df['count'].sum()
    _names = list(df.columns)
    _n = set(df[_names[0]])
    _color = alt.Color(f'{_names[0]}:N', title = _names[0].title(), legend=alt.Legend(orient=lgpos, title = None)).scale(scheme = 'turbo')
    return (
    alt.Chart(df)
    .transform_calculate('percent', f'datum.count / {total}')
    .mark_arc(innerRadius = 20)
    .encode(
        theta=alt.Theta('count:Q'),
        color = _color,
        tooltip = [
            alt.Tooltip(f'{_names[0]}:O', title = _names[0].title()),
            alt.Tooltip('percent:Q', title = "Percent", format = ".2%"),
            alt.Tooltip("count:Q", title = "Count")
        ]
    )
    .properties(
        height = 100, width = 150,
        title = alt.Title(title, align = 'center', anchor = 'middle')
    )
)

def convolutionbar(df: pd.DataFrame, title: str = ""):
    total = df['count'].sum()
    _names = list(df.columns)
    #_color = alt.Color(f'{_names[0]}:N', legend=None).scale(scheme = 'turbo')
    _color = alt.Color(f'{_names[0]}:N', legend=None).scale(scheme = 'lightgreyred')
    return (
        alt.Chart(df)
        .transform_calculate('percent', f'datum.count / {total}')
        .mark_bar()
        .encode(
            x=alt.X('count:Q', axis = None),
            y=alt.Y(f'{_names[0]}:N', title = None, sort = 'ascending'),
            color = _color,
            tooltip = [
                alt.Tooltip(f'{_names[0]}:O', title = _names[0].title()),
                alt.Tooltip('percent:Q', title = "Percent", format = ".2%"),
                alt.Tooltip("count:Q", title = "Count")
            ]
        )
        .properties(
            height = 100, width = 150,
            title = alt.Title(title, align = 'center', anchor = 'middle')
        )
)

def _makepanel(df, metric, title=None):
    '''Create an altair chart of read depth. Used internally by `depthplot()`'''
    nearest = alt.selection_point(
        nearest=True,
        on='mouseover',
        fields=['Position'],
        empty=False
    )

    base = (
        alt.Chart(df)
        .transform_fold(['Read Depth', 'Molecule Depth'], as_=['Metric', 'Depth'])
        .transform_filter(alt.datum.Metric == metric)
        .encode(
            x=alt.X('Position:Q', title='Position (Mb)')
                .scale(domainMin=0)
                .axis(labelExpr='datum.value / 1000000'),
            y=alt.Y('Depth:Q', title=None),
            color=alt.Color('Metric:N').legend(None),
        )
    )

    line = base.mark_line()

    points = base.mark_point(opacity=0).encode(
        tooltip=["Metric:N", 'Genomic Interval (bp):N', 'Depth:Q'],
        opacity=alt.condition(nearest, alt.value(1), alt.value(0))
    ).add_params(nearest)

    rule = base.mark_rule(color='gray', strokeWidth=1).encode(
        tooltip=["Metric:N", 'Genomic Interval (bp):N', 'Depth:Q'],
    ).transform_filter(nearest)

    return (
        alt.layer(line, points, rule)
        .properties(height=300, width=900, title=title or metric)
    )

def depthplot(df, title):
    alt.data_transformers.disable_max_rows()
    coverage = "Read Depth" in df.columns
    molcov = "Molecule Depth" in df.columns
    if coverage and molcov:
        res = alt.vconcat(
                _makepanel(df, 'Read Depth', title=f"{title} (Read Depth)"),
                _makepanel(df, 'Molecule Depth', title = f"{title} (Molecule Depth)")
            ).resolve_scale(x='independent', y='independent')
    else:
        _param = 'Read Depth' if coverage else 'Molecule Depth'
        res = _makepanel(df, _param, title = f"{title} ({_param})")
    return res.interactive(bind_y = False)

