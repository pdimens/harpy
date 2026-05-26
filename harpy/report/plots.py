from contextlib import contextmanager
import altair as alt
import pandas as pd
from harpy.report.theme import sv_colors

@contextmanager
def SafeRender(data):
    '''
    Given a Pandas/Polars DataFrame (or list), will temporarily set the Altair render engine to SVG
    to avoid the 5000 max row cutoff.
    '''
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
    alt.data_transformers.enable('default')
    return res.interactive(bind_y = False)

def sv_by_chromosome(variants: pd.DataFrame, title:str = ""):
    '''
    Return an Altair chart of SVs and their positions in a chromosome. This includes
    a dropdown selection to display a specific chromosome.
    '''
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
    '''
    Return an Altair chart of alignment depth in `window` bp intervals with a
    chromosome dropdown option that dynamically changes which chromosome's
    depths you see in the plot view.
    '''
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