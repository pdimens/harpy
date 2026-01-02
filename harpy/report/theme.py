import altair as alt

@alt.theme.register("harpy_theme", enable=True)
def harpy_theme() -> alt.theme.ThemeConfig:
    '''Altair theme settings for Harpy reports'''
    # Okabe Ito colorblindsafe palette
    # https://thenode.biologists.com/data-visualization-with-flying-colors/research/

    markColor = "#0072B2"  # blue, fifth color from category palette
    category_palette = ["#E69F00",
                        "#56B4E9",
                        "#009E73",
                        "#F0E442",
                        "#0072B2",
                        "#D55E00",
                        "#CC79A7",
                        "#000000"
                   ]
    sequential_palette = ["#cfe8f3", 
                          "#a2d4ec", 
                          "#73bfe2", 
                          "#46abdb", 
                          "#1696d2", 
                          "#12719e", 
                         ]
    return {
        "width": "container",
        "height": 400,
        "config": {
            "title": {
                "anchor": "start",
                "dy": -15,
                "fontSize": 18,
            },
            "arc": {"fill": markColor},
            "axis": {
                "domainColor": "#000000cc",
                "domainWidth": 0.5,
                "gridWidth": 0.2,
                "grid" : False,
                "titleFontWeight": 400,
                "titleFontSize": 15,
                "labelFontSize": 11,
                "labelColor": "#000000cc",
                "tickColor": "#000000cc",
                "tickWidth": 0.2,
                "titleColor": "#000000cc",
            },
            "axisBand": {"grid": False},
            "axisX": {"grid": True, "tickSize": 10},
            "axisY": {"domain": False, "grid": True, "tickSize": 0},
            "background": "#ffffff4d",
            "legend": {
                "labelFont": "sans-serif",
                "labelFontSize": 12,
                "symbolSize": 100,
                "symbolType": "square",
                "titleFontSize": 12,
                "titlePadding": 10,
                "title": "",
            },

           "range": {
               "category": category_palette,
               "diverging": sequential_palette,
           },
           "area": {
               "fill": markColor,
               "fillOpacity": 0.7,
               "line": True
           },
           "line": {
               "color": markColor,
               "stroke": markColor,
               "strokewidth": 5,
           },
           "trail": {
               "color": markColor,
               "stroke": markColor,
               "strokeWidth": 0,
               "size": 1,
           },
           "path": {
               "stroke": markColor,
               "strokeWidth": 0.5,
           },
           "point": {
               "filled": True,
           },
           "text": {
               "color": markColor,
               "fontSize": 11,
               "align": "right",
               "fontWeight": 400,
               "size": 11,
           }, 
           "bar": {
                "fill": markColor,
                "stroke": False,
            }, 
            "view": {
               "stroke": "none",
               "continuousWidth": 660,
               "continuousHeight": 400
           },
        }
    }