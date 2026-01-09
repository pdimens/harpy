import altair as alt
import pandas as pd

def by_chromosome_plot(variants: pd.DataFrame, title:str = ""):
    labels = variants['chromosome'].unique()
    input_dropdown = alt.binding_select(options=labels, name='Chromosome: ')
    selection = alt.selection_point("chrom_choice", fields=['chromosome'], value=labels[0], bind=input_dropdown)
    length_param = alt.param(expr='data("data_0")[0].length')
    highlight = alt.selection_point(name="highlight", on="pointerover", empty=False)
    stroke_color = (
        alt.when(highlight).then(alt.value("#7ae00d"))
        .otherwise(alt.value("transparent"))
    )

    return (
        alt.Chart(variants)
        .mark_bar(strokeWidth=1.5, cornerRadius=8)
        .encode(
            x=alt.X('start:Q')
                .scale(alt.Scale(domain=[0, length_param]))
                .axis(alt.Axis(title='Position (Mb)', labelExpr='datum.value / 1000000')),
            x2='end:Q',
            y=alt.Y('variant:N', title = "Variant Type"),
            color=alt.Color('variant:N').legend(None),
            tooltip=['variant:N', 'chromosome:N', 'start:Q', 'end:Q'],
            stroke=stroke_color
        )
        .transform_filter(selection)
        .add_params(selection, length_param, highlight)
        .properties(title= title)
    )

def depth_by_chromosome(records: pd.DataFrame, window: int = 50000, title:str = ""):
    '''
    Return an Altair chart of alignment depth in `window` bp intervals with a
    chromosome dropdown option that dynamically changes which chromosome's
    depths you see in the plot view.
    '''
    labels = records['contig'].unique()
    input_dropdown = alt.binding_select(options=labels, name='Contig: ')
    selection = alt.selection_point("chrom_choice", fields=['contig'], value=labels[0], bind=input_dropdown)
    length_param = alt.param(expr='max(pluck(data("data_0"), "position_end"))')
    highlight = alt.selection_point(name="highlight", on="pointerover", empty=False)
    stroke_color = (
        alt.when(highlight).then(alt.value("#7ae00d"))
        .otherwise(alt.value("transparent"))
    )

    return (
        alt.Chart(records)        
        .mark_bar(strokeWidth=2)
        .encode(
            x=alt.X('position:Q')
                .scale(alt.Scale(domain=[0, length_param]))
                .axis(alt.Axis(title='Position (Mb)', labelExpr='datum.value / 1000000')),
            y = 'count()',
            color = 'type:N',
            stroke = stroke_color
        )
        .transform_filter(selection)
        .add_params(selection, length_param, highlight)
        .properties(title= title)
        facet(row='key:N')
    )