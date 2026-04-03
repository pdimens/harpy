import altair as alt

# Okabe Ito colorblindsafe palette
# https://thenode.biologists.com/data-visualization-with-flying-colors/research/

markColor = "#0072B2"  # blue, fifth color from category palette
textColor = "#888888"

okabe_category = [
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
    "#000000"
]
okabe_sequential = [
    "#cfe8f3", 
    "#a2d4ec", 
    "#73bfe2", 
    "#46abdb", 
    "#1696d2", 
    "#12719e", 
]

def get_okabe(index:int|list[int], palette: str = "cat"):
    '''Return a color (or colors) from the `okabe_category` or `okabe_sequential` palettes by index'''
    if palette == "cat":
        if isinstance(index, int):
            return okabe_category[index]
        else:
            return [okabe_category[i] for i in index]
    else:
        if isinstance(index, int):
            return okabe_sequential[index]
        else:
            return [okabe_sequential[i] for i in index]


def palette(n_colors: int, palette: list[str]=[]):
    """
    Return a list of n equally-spaced colors from a `palette`. If `palette` is given, it must
    be a `list` of hex colors (e.g. `"#372466"`). If no palette is given, uses the
    256-color Turbo palette.

    Turbo palette from https://github.com/dbaranger/turbo_palette_R/blob/main/turbo_hex.txt
    """
    if not palette:
        # 256-color distribution of the Turbo palette
        palette = ["#30123B","#321543","#33184A","#341B51","#351E58","#36215F","#372466","#38276D","#392A73","#3A2D79","#3B2F80","#3C3286","#3D358B","#3E3891","#3F3B97","#3F3E9C","#4040A2","#4143A7","#4146AC","#4249B1","#424BB5","#434EBA","#4451BF","#4454C3","#4456C7","#4559CB","#455CCF","#455ED3","#4661D6","#4664DA","#4666DD","#4669E0","#466BE3","#476EE6","#4771E9","#4773EB","#4776EE","#4778F0","#477BF2","#467DF4","#4680F6","#4682F8","#4685FA","#4687FB","#458AFC","#458CFD","#448FFE","#4391FE","#4294FF","#4196FF","#4099FF","#3E9BFE","#3D9EFE","#3BA0FD","#3AA3FC","#38A5FB","#37A8FA","#35ABF8","#33ADF7","#31AFF5","#2FB2F4","#2EB4F2","#2CB7F0","#2AB9EE","#28BCEB","#27BEE9","#25C0E7","#23C3E4","#22C5E2","#20C7DF","#1FC9DD","#1ECBDA","#1CCDD8","#1BD0D5","#1AD2D2","#1AD4D0","#19D5CD","#18D7CA","#18D9C8","#18DBC5","#18DDC2","#18DEC0","#18E0BD","#19E2BB","#19E3B9","#1AE4B6","#1CE6B4","#1DE7B2","#1FE9AF","#20EAAC","#22EBAA","#25ECA7","#27EEA4","#2AEFA1","#2CF09E","#2FF19B","#32F298","#35F394","#38F491","#3CF58E","#3FF68A","#43F787","#46F884","#4AF880","#4EF97D","#52FA7A","#55FA76","#59FB73","#5DFC6F","#61FC6C","#65FD69","#69FD66","#6DFE62","#71FE5F","#75FE5C","#79FE59","#7DFF56","#80FF53","#84FF51","#88FF4E","#8BFF4B","#8FFF49","#92FF47","#96FE44","#99FE42","#9CFE40","#9FFD3F","#A1FD3D","#A4FC3C","#A7FC3A","#A9FB39","#ACFB38","#AFFA37","#B1F936","#B4F836","#B7F735","#B9F635","#BCF534","#BEF434","#C1F334","#C3F134","#C6F034","#C8EF34","#CBED34","#CDEC34","#D0EA34","#D2E935","#D4E735","#D7E535","#D9E436","#DBE236","#DDE037","#DFDF37","#E1DD37","#E3DB38","#E5D938","#E7D739","#E9D539","#EBD339","#ECD13A","#EECF3A","#EFCD3A","#F1CB3A","#F2C93A","#F4C73A","#F5C53A","#F6C33A","#F7C13A","#F8BE39","#F9BC39","#FABA39","#FBB838","#FBB637","#FCB336","#FCB136","#FDAE35","#FDAC34","#FEA933","#FEA732","#FEA431","#FEA130","#FE9E2F","#FE9B2D","#FE992C","#FE962B","#FE932A","#FE9029","#FD8D27","#FD8A26","#FC8725","#FC8423","#FB8122","#FB7E21","#FA7B1F","#F9781E","#F9751D","#F8721C","#F76F1A","#F66C19","#F56918","#F46617","#F36315","#F26014","#F15D13","#F05B12","#EF5811","#ED5510","#EC530F","#EB500E","#EA4E0D","#E84B0C","#E7490C","#E5470B","#E4450A","#E2430A","#E14109","#DF3F08","#DD3D08","#DC3B07","#DA3907","#D83706","#D63506","#D43305","#D23105","#D02F05","#CE2D04","#CC2B04","#CA2A04","#C82803","#C52603","#C32503","#C12302","#BE2102","#BC2002","#B91E02","#B71D02","#B41B01","#B21A01","#AF1801","#AC1701","#A91601","#A71401","#A41301","#A11201","#9E1001","#9B0F01","#980E01","#950D01","#920B01","#8E0A01","#8B0902","#880802","#850702","#810602","#7E0502","#7A0403"]
    
    if n_colors == 1:
        return [palette[0]]
    
    # Calculate step size to evenly space colors
    step = len(palette) / (n_colors - 1)
    
    # Select colors at regular intervals, ensuring last color is included
    selected_colors = [palette[min(int(i * step), len(palette) - 1)] for i in range(n_colors)]
    
    return selected_colors

def sv_colors(name) -> str:
    '''Return a color from that okabe palette that's mapped to a type of structural variant'''
    d = {
        "INV": "#56B4E9",
        "DEL": "#E69F00",
        "DUP": "#009E73",
        "BND": "#CC79A7"
    }
    return d.get(name, "#CC79A7")

@alt.theme.register("harpy_theme", enable=True)
def harpy_theme() -> alt.theme.ThemeConfig:
    '''Altair theme settings for Harpy reports'''

    return {
        "width": "container",
        "height": 400,
        "config": {
            "title": {
                "anchor": "start",
                "dy": -15,
                "fontSize": 18,
                "color": textColor,
                "subtitleColor": textColor,
                "subtitleFontSize": 14
            },
            "arc": {"fill": markColor},
            "axis": {
                "domainColor": textColor,
                "domainWidth": 0.5,
                "gridWidth": 0.2,
                "grid" : False,
                "titleFontWeight": 400,
                "titleFontSize": 15,
                "labelFontSize": 11,
                "labelColor": textColor,
                "tickColor": textColor,
                "tickWidth": 0.2,
                "titleColor": textColor,
            },
            "axisBand": {"grid": False},
            "axisX": {"grid": True, "tickSize": 10, "domainCap" : 'round'},
            "axisY": {"domain": False, "grid": True, "tickSize": 0, "domainCap" : 'round'},
            "background": "transparent",
            "header": {
                "labelFontSize": 15,
                "labelColor": textColor,
                "titleFontSize": 15,
                "titleFontWeight": 600,
                "titleColor": textColor
            },
            "legend": {
                "labelFont": "sans-serif",
                "labelFontSize": 12,
                "labelColor": textColor,
                "symbolSize": 100,
                "symbolType": "circle",
                "titleFontSize": 12,
                "titlePadding": 10,
                "titleColor": textColor,
                "title": None,
            },

           "range": {
               "category": okabe_category,
               "diverging": okabe_sequential,
           },
           "area": {
               "fill": markColor,
               "fillOpacity": 0.7,
               "line": True,
               "strokeOpacity": 0.7
           },
           "line": {
               "color": markColor,
               "stroke": markColor,
               "strokeWidth": 2,
           },
           "trail": {
               "color": markColor,
               "stroke": markColor,
               "strokeWidth": 0,
               "size": 1,
           },
           "point": {
               "filled": True,
           },
           "text": {
               "color": textColor,
               "fontSize": 11,
               "align": "right",
               "fontWeight": 400,
               "size": 11,
           }, 
           "bar": {
                "fill": markColor,
                "stroke": None,
            }, 
            "view": {
               "stroke": None,
               "continuousWidth": 660,
               "continuousHeight": 400
           },
        }
    }