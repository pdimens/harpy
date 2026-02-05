from datetime import datetime
from IPython.display import display, HTML
import itables
itables.options.warn_on_undocumented_option=False

class colored_boxes:
  '''
  The colored boxes class that draws boxes at the top of the reports.
  '''
  def __init__(self):
    self.boxes: list[str] = []

  def add(self, value, label, color = "#aeaeaeff", width: int = 200, height: int = 85, fontsize = 42):
    '''
    Return the html of a colored box object with `value` and `label`. Width given in pixels.
    '''
    if isinstance(value, str):
      _val = value
    else:
      _val = f"{value:,.2f}".rstrip('0').rstrip('.')
    _html = '''<div style="background-color: {}; width: {}; height: {}; display: flex;\
  flex-direction: column; align-items: center; justify-content: center;\
  color: white; border-radius: 15px; box-shadow: 0 4px 10px rgba(0,0,0,0.19);\
  padding: 15px; box-sizing: border-box;"><div style="font-size: 14px; font-weight: normal;\
    margin-bottom: 0px; margin-top: 3px; opacity: 0.9; text-transform: uppercase;\
  letter-spacing: 1px;">{}</div><div style="font-size: {}px; font-weight: bold;">{}</div></div>'''
    _html = _html.format(color, f"{width}px", f"{height}px", label, fontsize, _val)
    self.boxes.append(_html)
    return self

  def conditional(self, value, label, cutoff: int|float, lower_bad: bool = True, as_percent:bool = False, width:int = 200, height:int = 85):
    '''
    Return the html of a colored box object with `value` and `label`. Use `as_percent` to multiply
    the value by 100 for printing purposes.
    The `color` is either yellow or green depending on what is determined better or worse than the `cutoff`:
    - `lower_bad=True`: `color` = yellow when value < cutoff (default)
    - `lower_bad=False`: `color` = yellow when value >= cutoff 
    '''
    if lower_bad:
      color = "#f6ab3c" if value < cutoff else "#68ae6b"
    else:
      color = "#f6ab3c" if value >= cutoff else "#68ae6b"

    return self.add(value if not as_percent else f"{value * 100}%", label, color, width = width, height = height)

  def render(self, gap = 15):
    '''Display all the colored boxes stored in `self.boxes` in one continuous wrapped row'''
    html_content = '<div style="display: flex; flex-wrap: wrap; gap: {}px; width: 100%;">{}</div>'
    return display(HTML(html_content.format(gap, "\n".join(self.boxes))))

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

def standard_itable(data, filename:str, caption: str|None= None, fixedcols: int|None = None, coldefs = [], html: bool = False):
    '''
    Boilerplate version of `itables.show` to reduce in-report boilerplate. Use `filename` to specify the output prefix for
    exports, `caption` to set an optional caption, `fixcols` to freeze this-many columns on the left, and `coldefs` to
    populate `columnDefs` with wacky JS.
    '''
    return itables.show(
        data,
        caption = caption,
        buttons = [
            "pageLength", "colvis",
            {"extend": "csvHtml5", "title": filename},
            {"extend": "excelHtml5", "title": filename}
        ],
        fixedColumns = {"start": fixedcols} if fixedcols else {},
        ordering = {"indicators": False, "handler": False},
        columnControl = [["order", "searchDropdown"]],
        showIndex= False,
        scrollX =  True,
        autoWidth=True,
        searching = False,
        style = "width:100%",
        classes = "display nowrap compact",
        allow_html=html,
        columnDefs = coldefs
    )

def proportionbar(value, max_value):
    """
    Create an HTML bar for a given value. This is intended to create a
    proportion-bar formatting for values in a pandas column by way of
    ```
    max_val = df['ValidBarcodes'].max()

    df_styled['ValidBarcodes'] = df['ValidBarcodes'].apply(lambda x: proportionbar(x, max_val))
    ```
    """
    percentage = (value / max_value) * 100 if max_value > 0 else 0
    return f'''
        <div style="width: 100%; display: flex; align-items: center;">
            <div style="flex-grow: 1; background: linear-gradient(90deg, 
                #5cb85c {percentage}%, transparent {percentage}%); 
                padding: 2px 5px;">
                {value}
            </div>
        </div>
    '''
