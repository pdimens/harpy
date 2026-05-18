import base64
import html
import json
from pathlib import Path
from datetime import datetime
from io import BytesIO
from IPython.display import display, HTML, Image
from PIL import Image as PImage
import uuid

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

    def add(self, value, label, plus_minus: str|int|float = '', textcol=None, units=None):
        '''
        Return the html of a colored box object with `value` and `label`. Units get added to the plusminus text if they are present,
        otherwise added to the main value.
        
        Args:
            value: Main metric to display
            label: Title of the metric
            plus_minus: smaller text to add if there is a plus-minus value
            textcol: color to make the value text
            units: add units to the end without spaces (e.g. 'bp')
        '''
        _val = value if isinstance(value, str) else f"{value:,.2f}".rstrip('0').rstrip('.')
        _val = f"{_val}{units}" if units else _val
        _col = textcol or 'var(--tw-prose-body, #666666)'
        pm_str = f'±{plus_minus}' if plus_minus else ''
        self.boxes.append(
            self.html.format(label = label, value = _val, color = _col, plus_minus=pm_str)
        )
        return self

    def conditional(self, value, label, cutoff: int|float, lower_bad: bool = True, as_percent:bool = False, plus_minus: str|int|float = '', add_percent:bool = False, digits = None):
        '''
        Return the html of a colored box object with `value` and `label`. Use `as_percent` to multiply
        the value by 100 for printing purposes. Use `digits` as a stand-in for rounding. Using `add_percent` to add `%`
        the end without multiplying by 100.
        The `color` is either yellow or green depending on what is determined better or worse than the `cutoff`:
        - `lower_bad=True`: `color` = yellow when value < cutoff (default)
        - `lower_bad=False`: `color` = yellow when value >= cutoff
        '''
        if lower_bad:
            color = "#f6ab3c" if value < cutoff else None
        else:
            color = "#f6ab3c" if value >= cutoff else None

        if as_percent:
            _v = value * 100
            _v = round(_v, digits) if digits else _v
            value = f"{_v}%"
        elif add_percent:
            _v = round(value, digits) if digits else value
            value = f"{_v}%"
        if not isinstance(plus_minus, str):
            pm = round(plus_minus, digits) if digits else plus_minus
        else:
            pm = plus_minus
        return self.add(value, label, plus_minus=pm, textcol = color)

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
