from IPython.display import display, HTML
import itables

class colored_boxes:
  def __init__(self):
    '''
    The colored boxes class that draws boxes at the top of the reports.
    '''
    self.boxes = []

  def add(self, value, label, color = "#aeaeaeff", width = "200px", height = "90px"):
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

  def show(self, gap = 15):
    '''Display all the colored boxes stored in `self.boxes` in one continuous wrapped row'''
    html_content = '<div style="display: flex; flex-wrap: wrap; gap: {}px; width: 100%;">{}</div>'
    return display(HTML(html_content.format(gap, "\n".join(self.boxes))))

def display_multitext(*args):
    '''HTML-print all arguments with line breaks between them'''
    _html = "<p>{}<p>".format("<br>".join(str(i) for i in [*args]))
    return display(HTML(_html))

def standard_itable(data, filename:str, caption: str|None= None, fixedcols: int|None = None, coldefs = []):
    '''
    Boilerplate version of `itables.show` to reduce in-report boilerplate. Use `filename` to specify the output prefix for
    exports, `caption` to set an optional caption, `fixcols` to freeze this-many columns on the left, and `coldefs` to
    populate `columnDefs` with wacky JS.
    '''
    return itables.show(
        data,
        caption=caption,
        buttons=[
            "pageLength", "colvis",
            {"extend": "csvHtml5", "title": filename},
            {"extend": "excelHtml5", "title": filename}
        ],
        fixedColumns = {"start": fixedcols} if fixedcols else {},
        ordering={"indicators": False, "handler": False},
        columnControl=[["order", "searchDropdown"]],
        showIndex=False,
        scrollX = True,
        columnDefs=coldefs
    )
    