import altair as alt
import pandas as pd

def by_chromosome_plot(variants: pd.DataFrame, title:str = ""):
    labels = variants['chromosome'].unique()
    input_dropdown = alt.binding_select(options=labels, name='Chromosome ')
    selection = alt.selection_point(fields=['chromosome'], value=labels[0], bind=input_dropdown)
    length_param = alt.param(expr='data("data_0")[0].length')
    return (
        alt.Chart(variants)
        .mark_rect(strokeWidth=2, cornerRadius=8)
        .encode(
            x=alt.X('start:Q')
                .scale(alt.Scale(domain=[0, length_param])),
            x2='end:Q',
            y='variant:N',
            color='variant:N',
            tooltip=['variant:N', 'chromosome:N', 'start:Q', 'end:Q']
        )
        .add_params(selection)
        .transform_filter(selection)
        .add_params(length_param)
        .properties(title=title)
    )