import altair as alt
import numpy as np
import pandas as pd
from harpy.report.theme import sv_colors

def radian_positions(chr_lens: dict[str,int]) -> pd.DataFrame:
    '''
    Convert the dict of chromosome names and lengths to a DataFrame
    that's in polar (radian) coordinates.
    '''
    # Create dataframe
    chrom_df = pd.DataFrame(list(chr_lens.items()), columns=['chromosome', 'length'])

    # Calculate cumulative positions (where each chromosome starts)
    chrom_df['cum_start'] = np.insert(np.cumsum(chrom_df['length'].values)[:-1], 0, 0)
    chrom_df['cum_end'] = chrom_df['cum_start'] + chrom_df['length']

    # Total genome length
    total_length = chrom_df['length'].sum()

    # Convert to radians, 
    #scale_factor = 0.95 # leaving a 5% gap for visual separation
    chrom_df['theta_start'] = (chrom_df['cum_start'] / total_length) * 2 * np.pi
    chrom_df['theta_end'] = (chrom_df['cum_end'] / total_length) * 2 * np.pi
    return chrom_df

def circos_grid(chr_df: pd.DataFrame, rings: int):
    '''
    Given the polar-converted chromosome grid `chr_df` (created by `radian_positions`),
    return the circos-style base plot with chromosome names and grides.
    '''
    gap = 4
    radius_width = 35
    outer_radius = 240
    inner_radius = outer_radius - radius_width
    # Create the circular chromosome ideogram
    base_chart = (
        alt.Chart(chr_df, width = 650, height = 650)
        .mark_arc(
            outerRadius= outer_radius,
            innerRadius= inner_radius,
            stroke='#bebebe',
            color='#f4f4f4',
            filled=True,
            fillOpacity=1,
            strokeJoin="round",
            strokeWidth=1,
            padAngle=0.01
        )
        .encode(
            theta=alt.Theta('theta_start:Q')
                .scale(alt.Scale(domain=[0, 2*np.pi])),
            theta2='theta_end:Q'
        )
    )
    for _ in range(rings-1):
        outer_radius = inner_radius - gap
        inner_radius = outer_radius - (radius_width * 0.9)
        base_chart += (
        alt.Chart(chr_df)
        .mark_arc(
            outerRadius=outer_radius,
            innerRadius=inner_radius,
            stroke='#bebebe',
            color = '#f4f4f4',
            strokeJoin="round",
            filled=True,
            fillOpacity=1,
            strokeWidth=1,
            padAngle=0.01
        )
        .encode(
            theta=alt.Theta('theta_start:Q')
                .scale(alt.Scale(domain=[0, 2*np.pi])),
            theta2='theta_end:Q'
        )
    )

    # Add chromosome labels
    labels = (
        alt.Chart(chr_df)
        .transform_calculate(theta_mid='(datum.theta_start + datum.theta_end) / 2')
        .mark_text(radius=270, size=14, fill = "black")
        .encode(
            theta='theta_mid:Q',
            text='chromosome:N'
        )
    )

    return (
        base_chart + labels
        .properties(title='Structural Variants')
    )

def variant_bands(var_df: pd.DataFrame, baseplot):
    gap = 4
    radius_width = 35
    outer_radius = 240
    inner_radius = outer_radius - radius_width
    var_groups = var_df.groupby('variant')
    for variant, gdf in var_groups:
        _color = sv_colors(variant)
        baseplot += (
            alt.Chart(gdf)
            .mark_arc(
                outerRadius=outer_radius - 3,
                innerRadius=inner_radius + 3,
                stroke=_color,
                strokeWidth=2,
                color = _color,
                thetaOffset=0.01,
                theta2Offset=0.01
            )
        .encode(
            theta=alt.Theta('theta_start:Q', scale=alt.Scale(domain=[0, 2*np.pi])),
            theta2='theta_end:Q',
            tooltip=['variant:N', 'chromosome:N', 'start:Q', 'end:Q']
        )
            )
        outer_radius = inner_radius - gap
        inner_radius = outer_radius - (radius_width * 0.9)

    return baseplot

def plot_circos_sv(chr_lengths: dict[str,int], sv: pd.DataFrame):
    '''
    The all-in-one function that takes a dict of chromosome names:lengths
    and a DataFrame of structural variants and does all the transformations
    to return a circos-equivalent plot using Altair:
    | chromosome | start | end | variant |
    |:------------|-------|-----|---------|
    | chr1 | 1000000 | 2000000 | INV |
    | chr1 | 5000000 | 6500000 | DEL |
    '''
    chrom_df = radian_positions(chr_lengths)
    total_length = chrom_df['length'].sum()

    # Process inversions: merge with chromosome cumulative positions
    var_df = sv.merge(chrom_df[['chromosome', 'cum_start']], on='chromosome')

    # Calculate absolute genomic positions
    var_df['abs_start'] = var_df['cum_start'] + var_df['start']
    var_df['abs_end'] = var_df['cum_start'] + var_df['end']

    # Convert to radians
    var_df['theta_start'] = (var_df['abs_start'] / total_length) * 2 * np.pi
    var_df['theta_end'] = (var_df['abs_end'] / total_length) * 2 * np.pi

    # create plot base/grid
    plot = circos_grid(chrom_df, var_df['variant'].nunique())
    return (plot + variant_bands(var_df, plot)).configure(background = "white")