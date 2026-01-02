from typing import Any, Self


from datetime import datetime
from IPython.display import display, HTML
import itables

class colored_boxes:
  '''
  The colored boxes class that draws boxes at the top of the reports.
  '''
  def __init__(self):
    self.boxes: list[str] = []

  def add(self, value, label, color = "#aeaeaeff", width = "200px", height = "90px") -> Self:
    '''
    Return the html of a colored box object with `value` and `label`
    '''
    _val = f"{value:,.2f}".rstrip('0').rstrip('.')
    _html = '''<div style="background-color: {}; width: {}; height: {}; display: flex;\
  flex-direction: column; align-items: center; justify-content: center;\
  color: white; border-radius: 10px; box-shadow: 0 4px 12px rgba(0,0,0,0.15);\
  padding: 20px; box-sizing: border-box;"><div style="font-size: 14px; font-weight: normal;\
    margin-bottom: 2px; margin-top: 13px; opacity: 0.9; text-transform: uppercase;\
  letter-spacing: 1px;">{}</div><div style="font-size: 43px; font-weight: bold;">{}</div></div>'''
    _html = _html.format(color, width, height, label, _val)
    self.boxes.append(_html)
    return self

  def render(self, gap = 15):
    '''Display all the colored boxes stored in `self.boxes` in one continuous wrapped row'''
    html_content = '<div style="display: flex; flex-wrap: wrap; gap: {}px; width: 100%;">{}</div>'
    return display(HTML(html_content.format(gap, "\n".join(self.boxes))))

def print_time(*args):
    '''HTML-print the time and all arguments with line breaks between them'''
    _now = datetime.now().strftime("🗓️ %d %B, %Y 🕔 %H:%M")
    _html = "<p>{}<p>".format("<br>".join([_now] + [str(i) for i in [*args]]))
    return display(HTML(_html))

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
