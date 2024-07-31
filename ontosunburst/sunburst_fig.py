import os.path
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from ontosunburst.data_table_tree import *

# ==================================================================================================
# CONSTANTS
# ==================================================================================================

# Kwargs
C_MIN = 'c_min'
C_MAX = 'c_max'
C_MID = 'c_mid'
MAX_DEPTH = 'max_depth'
COLORSCALE = 'colorscale'
TITLE = 'title'
COLORBAR_LEGEND = 'colorbar_legend'
BG_COLOR = 'bg_color'
FONT_COLOR = 'font_color'
FONT_SIZE = 'font_size'
TABLE_TITLE = 'table_title'
TABLE_LEGEND = 'table_legend'
TABLE_COLOR = 'table_color'


# ==================================================================================================
# FUNCTIONS
# ==================================================================================================

# Figure creation
# --------------------------------------------------------------------------------------------------
def get_fig_kwargs(output: str, analysis: str, **kwargs):
    """ Generate a Sunburst figure and save it to output path.

        Parameters
        ----------
        output: str (optional, default=None)
            Path to output to save the figure without extension
        analysis: str (optional, default=topology)
            Analysis mode : topology or enrichment
        """
    def_colorscale = {TOPOLOGY_A: 'Viridis',
                      ENRICHMENT_A: 'RdBu'}
    def_titles = {TOPOLOGY_A: f'{os.path.basename(output)} : Proportion of classes',
                  ENRICHMENT_A: f'{os.path.basename(output)} : Classes enrichment representation'}
    def_colorbar = {TOPOLOGY_A: 'Count',
                    ENRICHMENT_A: 'Log10(p-value)'}
    def_c_min = {TOPOLOGY_A: 1, ENRICHMENT_A: -10}
    def_c_max = {TOPOLOGY_A: None, ENRICHMENT_A: 10}
    def_c_mid = {TOPOLOGY_A: None, ENRICHMENT_A: 0}

    c_min = kwargs.get(C_MIN, def_c_min[analysis])
    c_max = kwargs.get(C_MAX, def_c_max[analysis])
    c_mid = kwargs.get(C_MID, def_c_mid[analysis])
    max_depth = kwargs.get(MAX_DEPTH, 7)
    colorscale = px.colors.get_colorscale(kwargs.get(COLORSCALE, def_colorscale[analysis]))
    title = kwargs.get(TITLE, def_titles[analysis])
    colorbar_legend = kwargs.get(COLORBAR_LEGEND, def_colorbar[analysis])
    background_color = kwargs.get(BG_COLOR, 'rgba(255, 255, 255, 0)')
    font_color = kwargs.get(FONT_COLOR, '#111111')
    font_size = kwargs.get(FONT_SIZE, 20)
    table_title = kwargs.get(TABLE_TITLE, 'Significant p-values')
    table_legend = kwargs.get(TABLE_LEGEND, 'IDs')
    table_color = kwargs.get(TABLE_COLOR, '#666666')

    return c_min, c_max, c_mid, max_depth, colorscale, title, colorbar_legend, background_color, \
        font_color, font_size, table_title, table_legend, table_color


def generate_sunburst_fig(data: DataTable, output: str, analysis: str = TOPOLOGY_A,
                          test=BINOMIAL_TEST, significant: Dict = None, ref_set: bool = True,
                          write_fig: bool = True, **kwargs) -> go.Figure:
    """ Generate a Sunburst figure and save it to output path.

    Parameters
    ----------
    data: DataTable

    output: str
        Path to output to save the figure without extension
    analysis: str (optional, default=topology)
        Analysis mode : topology or enrichment
    test: str (optional, default=Binomial)
        Type of test for enrichment analysis : Binomial or Hypergeometric
    root_cut: str (optional, default=ROOT_CUT)
        mode for root cutting (uncut, cut, total)
    ref_set: bool (optional, default=True)
    write_fig: bool (optional, default=True)
        True to write the html figure, False to only return figure
    **kwargs

    Returns
    -------
    go.Figure
    """
    c_min, c_max, c_mid, max_depth, colorscale, title, colorbar_legend, background_color, \
        font_color, font_size, table_title, table_legend, table_color = \
        get_fig_kwargs(output, analysis, **kwargs)

    if analysis == TOPOLOGY_A:
        fig = go.Figure(go.Sunburst(labels=data.labels, parents=data.parents,
                                    values=data.relative_prop, ids=data.ids,
                                    hoverinfo='label+text', maxdepth=max_depth,
                                    branchvalues='total',
                                    hovertext=get_hover_fig_text(data, TOPOLOGY_A, ref_set),
                                    marker=dict(colors=data.count, colorscale=colorscale,
                                                cmin=c_min, cmax=c_max, cmid=c_mid, showscale=True,
                                                colorbar=dict(title=dict(text=colorbar_legend)))))
        fig.update_layout(title=dict(text=title, x=0.5, xanchor='center'))

    elif analysis == ENRICHMENT_A:
        fig = make_subplots(rows=1, cols=2,
                            column_widths=[0.3, 0.7],
                            vertical_spacing=0.03,
                            subplot_titles=(table_title, title),
                            specs=[[{'type': 'table'}, {'type': 'sunburst'}]])

        fig.add_trace(go.Sunburst(labels=data.labels, parents=data.parents,
                                  values=data.relative_prop, ids=data.ids,
                                  hovertext=get_hover_fig_text(data, ENRICHMENT_A, ref_set),
                                  hoverinfo='label+text', maxdepth=max_depth,
                                  branchvalues='total',
                                  marker=dict(colors=data.p_val, colorscale=colorscale,
                                              cmid=c_mid, cmax=c_max, cmin=c_min, showscale=True,
                                              colorbar=dict(title=dict(text=colorbar_legend)))),
                      row=1, col=2)

        fig.add_trace(go.Table(header=dict(values=[table_legend, f'{test} test P-value'],
                                           fill=dict(color=table_color), height=40,
                                           font=dict(size=font_size)),
                               cells=dict(values=[list(significant.keys()),
                                                  list(significant.values())],
                                          fill=dict(color=table_color), height=35,
                                          font=dict(size=font_size*0.80))),
                      row=1, col=1)
    else:
        raise ValueError('Wrong type input')
    fig.update_layout(paper_bgcolor=background_color, font_color=font_color, font_size=font_size)
    fig.update_annotations(font_size=font_size*1.5)
    if write_fig:
        fig.write_html(f'{output}.html')
    return fig


def get_hover_fig_text(data: DataTable, analysis: str, ref_set: bool) \
        -> List[str]:
    """

    Parameters
    ----------
    data
    analysis
    ref_set

    Returns
    -------

    """
    if analysis == ENRICHMENT_A:
        return [f'P value: {10 ** (-data.p_val[i])}<br>'
                f'{COUNT}: <b>{data.count[i]}</b><br>'
                f'{REF_COUNT}: {data.ref_count[i]}<br>'
                f'{PROP}: <b>{round(data.prop[i] * 100, 2)}%</b><br>'
                f'{REF_PROP}: {round(data.ref_prop[i] * 100, 2)}%<br>'
                f'{IDS}: {data.onto_ids[i]}'
                if data.p_val[i] > 0 else
                f'P value: {10 ** data.p_val[i]}<br>'
                f'{COUNT}: <b>{data.count[i]}</b><br>'
                f'{REF_COUNT}: {data.ref_count[i]}<br>'
                f'{PROP}: <b>{round(data.prop[i] * 100, 2)}%</b><br>'
                f'{REF_PROP}: {round(data.ref_prop[i] * 100, 2)}%<br>'
                f'{IDS}: {data.onto_ids[i]}'
                for i in range(data.len)]
    elif analysis == TOPOLOGY_A:
        if ref_set:
            return [f'{COUNT}: <b>{data.count[i]}</b><br>'
                    f'{REF_COUNT}: {data.ref_count[i]}<br>'
                    f'{PROP}: <b>{round(data.prop[i] * 100, 2)}%</b><br>'
                    f'{REF_PROP}: {round(data.ref_prop[i] * 100, 2)}%<br>'
                    f'{IDS}: {data.onto_ids[i]}'
                    for i in range(data.len)]
        else:
            return [f'{COUNT}: <b>{data.count[i]}</b><br>'
                    f'{PROP}: <b>{round(data.prop[i] * 100, 2)}%</b><br>'
                    f'{IDS}: {data.onto_ids[i]}'
                    for i in range(data.len)]
