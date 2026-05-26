from contextlib import contextmanager
import altair as alt
import pandas as pd
from harpy.report.theme import sv_colors

@contextmanager
def SafeRender(data):
    """
    Context manager that temporarily switches Altair to SVG rendering and disables the 5000-row cutoff when given large data.
    
    When entered, saves Altair's current renderer. If len(data) > 5000, disables Altair's max-rows transformer and enables the "svg" renderer. On exit, restores the default data transformer with max_rows=5000 and re-enables the previously active renderer.
    
    Parameters:
        data (Sequence | pandas.DataFrame | polars.DataFrame): Object with a defined length used to decide whether to switch renderers; if its length is greater than 5000 the renderer and transformer are adjusted.
    """
    previous = alt.renderers.active
    if len(data) > 5000:
        alt.data_transformers.disable_max_rows()
        alt.renderers.enable("svg")
    try:
        yield
    finally:
        alt.data_transformers.enable('default', max_rows=5000)
        alt.renderers.enable(previous)

def piechart(df: pd.DataFrame, title: str = "", goodval: str|None = None, lgtitle: bool = False, lgpos: str = "bottom"):
    """
    Create a small Altair pie chart from a dataframe whose first column is the categorical field and which contains a 'count' column.
    
    Parameters:
        df (pd.DataFrame): DataFrame where the first column is the categorical field and 'count' holds numeric counts.
        title (str): Chart title shown centered above the chart.
        goodval (str | None): If provided, that category is colored green and the first other category is colored gray; otherwise colors are assigned automatically.
        lgtitle (bool): If True the legend displays a title (the categorical field name); if False the legend title is removed.
        lgpos (str): Legend orientation (e.g., 'bottom', 'right').
    
    Returns:
        alt.Chart: An Altair pie chart encoding category color, arc theta by 'count', and tooltips for category, percent (formatted as `.2%`), and count; sized to 150×100 with a centered title.
    """
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
    """
    Builds a small donut-style pie chart showing category proportions from a count DataFrame.
    
    Parameters:
        df (pd.DataFrame): Input table containing a categorical column (first column) and a numeric 'count' column; proportions are computed from 'count'.
        title (str): Chart title; shown centered.
        lgpos (str): Legend orientation (e.g., 'bottom', 'right').
    
    Returns:
        alt.Chart: An Altair arc chart (innerRadius=20) encoding slice size by `count`, color by the first DataFrame column, and tooltips for category, percent, and count.
    """
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
    """
    Builds a compact horizontal bar chart of categorical counts with percent and count tooltips.
    
    Parameters:
        df (pandas.DataFrame): DataFrame that must contain a 'count' column; the first DataFrame column is used as the categorical field.
        title (str): Chart title displayed centered; defaults to an empty string.
    
    Returns:
        alt.Chart: An Altair bar chart sized 150x100 that encodes categories on the y-axis, count on the x-axis, colors bars by category (lightgreyred scheme), and provides tooltips for category, percent (as `.2%`), and count.
    """
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
    """
    Create an interactive Altair panel showing depth for a single metric.
    
    Parameters:
        df (pandas.DataFrame): Source dataframe containing at least 'Position', 'Read Depth', and 'Molecule Depth' columns.
        metric (str): Name of the metric to display; expected values include 'Read Depth' or 'Molecule Depth'.
        title (str | None): Optional chart title; if omitted the metric name is used.
    
    Returns:
        alt.Chart: A layered Altair chart combining a line, hoverable points, and a hover rule, sized to 900×300.
    """
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
    """
    Builds an interactive depth chart showing read depth and/or molecule depth.
    
    If the dataframe contains both "Read Depth" and "Molecule Depth" columns, the chart vertically stacks two panels (one per metric) with independent x and y scales; otherwise it renders a single panel for the available metric. The resulting chart is interactive with vertical (y) panning/zooming disabled.
    
    Parameters:
        df (pandas.DataFrame): Input dataframe that must contain one or both of the columns "Read Depth" and "Molecule Depth".
        title (str): Base title used for the chart and panel subtitles.
    
    Returns:
        alt.Chart: An interactive Altair chart displaying depth information; if both metrics are present, a vertically concatenated chart of two panels is returned, otherwise a single-panel chart is returned.
    """
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

def sv_by_chromosome(variants: pd.DataFrame, title:str = ""):
    """
    Render an interactive Altair chart of structural variants filtered by chromosome.
    
    Parameters:
        variants (pd.DataFrame): DataFrame of structural variants. Must contain columns
            'Contig', 'Start', 'End', 'Type', and optionally 'N Samples' and 'Samples'.
        title (str): Chart title prefix; if empty the chart uses a dynamic title based on the selected contig.
    
    Returns:
        alt.Chart: An Altair chart that displays variant spans across genomic positions with a contig dropdown,
            hover highlighting, zoomable x-axis, and tooltips for type, contig, start, end, length, and sample counts.
    """
    _sv  = ["Inversion", "Deletion", "Duplication", "Breakend"]
    _col = [sv_colors(i) for i in ['INV', 'DEL', 'DUP', 'BND']]

    labels = variants['Contig'].unique()
    input_dropdown = alt.binding_select(options=labels, name='Contig: ')
    selection = alt.selection_point(name = "chrom_choice", fields=['Contig'], value=labels[0], bind=input_dropdown)
    length_param = alt.param(expr='data("data_0")[0].length')
    highlight = alt.selection_point(name="highlight", on="pointerover", empty=False)
    zoom = alt.selection_interval(bind='scales', encodings=['x'])
    stroke_color = (
        alt.when(highlight)
        .then(alt.value("#7ae00d"))
        .otherwise(alt.Color('Type:N').scale(domain = _sv, range = _col))
    )
    dynamic_title = alt.Title(alt.expr(f'"Structural Variants on " + {selection.name}.Contig'), subtitle = "Variants should be considered putative")

    return (
        alt.Chart(variants)
        .transform_calculate(var_length = 'datum.End - datum.Start')
        .transform_filter(selection)
        .mark_bar(strokeWidth = 2, cornerRadius=8, opacity = 0.75)
        .encode(
            x=alt.X('Start:Q')
                .scale(domain=[0, length_param])
                .axis(title='Position (Mb)', labelExpr='datum.value / 1000000'),
            x2='End:Q',
            y=alt.Y('Type:N', title = "Variant Type"),
            color=alt.Color('Type:N', legend = None)
                .scale(domain = _sv, range = _col),
            tooltip=[
                alt.Tooltip('Type:N', title = "Variant Type"),
                alt.Tooltip('Contig:N', title = "Contig"),
                alt.Tooltip('Start:Q', title = "Start", format = ','),
                alt.Tooltip('End:Q', title = "End", format = ','),
                alt.Tooltip('var_length:Q', title = "Length", format = ','),
                alt.Tooltip('N Samples:Q', title = "# Samples"),
                alt.Tooltip('Samples:N', title = "Samples")
            ],
            stroke=stroke_color
        )
        .add_params(selection, length_param, highlight, zoom)
        .properties(title= dynamic_title)
    )

def depth_by_chromosome(records: pd.DataFrame, title:str = ""):
    """
    Create an interactive Altair chart showing alignment depth by position with a contig dropdown filter.
    
    The chart is faceted by the `key` column, uses a dropdown bound to the distinct values of `records['contig']` to select which contig to display, and scales the x-axis domain to the maximum `position_end` found in the plotted data. Bars are colored by `type` and receive a green stroke on pointerover.
    
    Parameters:
        records (pd.DataFrame): DataFrame containing at minimum the columns `contig`, `position`, `position_end`, `type`, and `key`.
        title (str): Optional chart title.
    
    Returns:
        alt.Chart: An Altair chart faceted by `key` with a contig dropdown, hover highlighting, and x-axis labeled in megabases.
    """
    labels = records['contig'].unique()
    input_dropdown = alt.binding_select(options=labels, name='Contig: ')
    selection = alt.selection_point(name = "chrom_choice", fields=['contig'], value=labels[0], bind=input_dropdown)
    length_param = alt.param(expr='max(pluck(data("data_0"), "position_end"))')
    highlight = alt.selection_point(name="highlight", on="pointerover", empty=False)
    stroke_color = (
        alt.when(highlight)
        .then(alt.value("#7ae00d"))
        .otherwise(alt.value("transparent"))
    )
    return (
        alt.Chart(records)        
        .mark_bar(strokeWidth=2)
        .encode(
            x=alt.X('position:Q')
                .scale(domain=[0, length_param])
                .axis(title='Position (Mb)', labelExpr='datum.value / 1000000'),
            y = 'count()',
            color = 'type:N',
            stroke = stroke_color
        )
        .transform_filter(selection)
        .add_params(selection, length_param, highlight)
        .properties(title= title)
        .facet(row='key:N')
    )