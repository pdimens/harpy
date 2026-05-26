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
        """
        Initialize a StatsBox builder with sizing and style defaults and prepare an HTML template.
        
        Parameters:
            height (int): Pixel height for each stat box (used in the template).
            label_fontsize (int): Font size in pixels for the label text.
            value_fontsize (int): Font size in pixels for the main value text.
            text_gap (int): Vertical gap in pixels between the label and value.
        
        Notes:
            - Creates an empty list `self.boxes` to accumulate rendered box HTML snippets.
            - Computes `self.pm_fontsize` as one third of `value_fontsize`.
            - Builds `self.html`, a reusable HTML template with placeholders `{label}`, `{value}`, `{plus_minus}`, and `{color}`.
        """
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
        """
        Add a stat box entry with a label, formatted value, optional plus/minus text, and optional text color.
        
        Parameters:
            value (str|int|float): Main metric to display. Non-string numeric values are formatted with two decimals and trailing zeros/decimal point removed.
            label (str): Title for the metric shown above the value.
            plus_minus (str|int|float, optional): Supplemental smaller text displayed prefixed with '±' when provided.
            textcol (str, optional): CSS color value to use for the value text; defaults to the module's prose color when omitted.
            units (str, optional): Unit suffix appended directly to the formatted value when provided (no added space).
        
        Returns:
            self: The same StatsBox instance to allow method chaining.
        """
        _val = value if isinstance(value, str) else f"{value:,.2f}".rstrip('0').rstrip('.')
        _val = f"{_val}{units}" if units else _val
        _col = textcol or 'var(--tw-prose-body, #666666)'
        pm_str = f'±{plus_minus}' if plus_minus else ''
        self.boxes.append(
            self.html.format(label = label, value = _val, color = _col, plus_minus=pm_str)
        )
        return self

    def conditional(self, value, label, cutoff: int|float, lower_bad: bool = True, as_percent:bool = False, plus_minus: str|int|float = '', add_percent:bool = False, digits = None):
        """
        Add a stat box for the given value and label, applying a conditional text color based on a cutoff.
        
        Parameters:
            value: The numeric or string value to display. If numeric and formatted as percent (see below), it will be converted to a string with a trailing '%' as appropriate.
            label: The label text for the stat box.
            cutoff: Threshold used to decide conditional coloring.
            lower_bad (bool): If True, a value lower than `cutoff` is considered bad; if False, a value greater than or equal to `cutoff` is considered bad.
            as_percent (bool): If True, multiply a numeric `value` by 100, optionally round using `digits`, and append '%'.
            plus_minus (str|int|float): Optional secondary value displayed alongside `value`. If numeric and `digits` is provided, it will be rounded.
            add_percent (bool): If True and `as_percent` is False, append '%' to the (optionally rounded) numeric `value` without multiplying by 100.
            digits (int|None): Number of decimal places to round numeric `value` and `plus_minus` before formatting; if None, no rounding is applied.
        
        Returns:
            self: The same StatsBox instance (allows method chaining).
        
        Behavior:
            Sets the text color to yellow (#f6ab3c) when the value is considered bad according to `lower_bad` and `cutoff`:
            - If `lower_bad` is True: yellow when value < cutoff.
            - If `lower_bad` is False: yellow when value >= cutoff.
            When not yellow, the default text color is used.
        """
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
        """
        Render accumulated stat boxes as a horizontal flex row and display them in the notebook.
        
        Parameters:
            gap (int): Pixel gap between boxes; boxes will wrap to multiple rows if they exceed the container width (default: 5).
        """
        container_html = f'''<div style="display: flex; gap: {gap}px; 
            flex-wrap: wrap;">{" ".join(self.boxes)}</div>'''
        display(HTML(container_html))

def print_time(*args):
    """
    Display the current date and time followed by provided arguments as an HTML paragraph with line breaks.
    
    Parameters:
        *args: Any values to display after the timestamp; each value is converted to a string and separated by an HTML line break (`<br>`).
    
    Notes:
        The timestamp is formatted as "🗓️ DD Month, YYYY 🕔 HH:MM".
    """
    _now = datetime.now().strftime("🗓️ %d %B, %Y 🕔 %H:%M")
    _html = "<p>{}<p>".format("<br>".join([_now] + [str(i) for i in [*args]]))
    return display(HTML(_html))

def print_html(*args):
    """
    Display provided values as a single HTML paragraph, with each value converted to a string and separated by HTML line breaks.
    
    Parameters:
        *args: Values to render; each value is stringified and joined with `<br>` between items.
    """
    display(HTML(
        "<p>{}</p>".format("<br>".join(str(i) for i in [*args]))
    ))

def embed_image(x: str, scale: float = 0.5):
    """
    Load an image file, optionally downscale it, and return an IPython.display.Image with the image embedded.
    
    Parameters:
        x (str): Path to the image file.
        scale (float): Scaling multiplier; values less than 1 will downscale the image, values >= 1 leave the image size unchanged.
    
    Returns:
        IPython.display.Image: An Image object containing the encoded image data. The image is encoded as PNG if the filename contains "png" (case-insensitive), otherwise as JPEG.
    """
    image = PImage.open(x)
    if scale < 1:  
        image = image.resize((int(image.size[0] * scale), int(image.size[1] * scale)), PImage.LANCZOS)  
    buf = BytesIO()  
    image.save(buf, format = "PNG" if "png" in x.lower() else "jpeg")  
    return Image(data=buf.getvalue(), format="png" if "png" in x.lower() else "jpeg", embed=True)

def image_viewer(label: str, image_dir: str, pattern: str, sortkey = None, thing_to_select: str = "sample", recursive: bool = False, scale: float = 1.0, option_key = None):  
    """
    Display an HTML/JavaScript image viewer in a Jupyter notebook that embeds matching image files as base64 data URIs and provides a dropdown to select which image to show.
    
    Parameters:
        label (str): Text label shown next to the dropdown selector.
        image_dir (str): Directory to search for image files.
        pattern (str): Glob pattern to match image filenames (e.g., "*.png" or "**/*.jpg").
        sortkey (callable | None): Optional key function passed to sorted() to control option order.
        thing_to_select (str): Noun used in the selector placeholder (e.g., "sample"); purely for display.
        recursive (bool): If True, search files recursively (rglob); otherwise use non-recursive glob.
        scale (float): If less than 1.0, downscale images by this factor before embedding.
        option_key (callable | None): Function mapping a Path to the option name shown and used as the selector value.
            Defaults to p.stem for non-recursive mode or p.parents[1].name for recursive mode.
    """
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
