import altair as alt
import pandas as pd
from harpy.report.theme import sv_colors

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
        .mark_bar(strokeWidth = 2, cornerRadius=8, opacity = 0.7)
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