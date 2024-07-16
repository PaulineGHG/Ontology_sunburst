import os.path
from typing import List, Dict
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


def generate_sunburst_fig(data: Dict[str, List[str or int or float]], output: str,
                          analysis: str = TOPOLOGY_A, ref_classes_abundance=None,
                          test=BINOMIAL_TEST, total: bool = True,
                          root_cut: str = ROOT_CUT, ref_base: bool = True, **kwargs) -> go.Figure:
    """ Generate a Sunburst figure and save it to output path.

    Parameters
    ----------
    data: Dict[str, List[str or int or float]]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
            - reference abundance value : Reference_count (int)
            - proportion : Proportion (0 < float <= 1)
            - reference proportion : Reference_proportion (0 < float <= 1)
            - branch proportion : Relative_prop
    output: str (optional, default=None)
        Path to output to save the figure without extension
    analysis: str (optional, default=topology)
        Analysis mode : topology or enrichment
    ref_classes_abundance: Dict[str, int] (optional, default=None)
        Abundances of reference set classes
    test: str (optional, default=Binomial)
        Type of test for enrichment analysis : Binomial or Hypergeometric
    total: bool (optional, default=True)
        True to have branch values proportional of the total parent
    root_cut: str (optional, default=ROOT_CUT)
        mode for root cutting (uncut, cut, total)
    ref_base: bool (optional, default=True)
    **kwargs

    Returns
    -------
    go.Figure
    """
    c_min, c_max, c_mid, max_depth, colorscale, title, colorbar_legend, background_color, \
        font_color, font_size, table_title, table_legend, table_color = \
        get_fig_kwargs(output, analysis, **kwargs)

    if total:
        branch_values = 'total'
        values = data[RELAT_PROP]
    else:
        branch_values = 'remainder'
        values = data[COUNT]

    if analysis == TOPOLOGY_A:
        data = data_cut_root(data, root_cut)
        fig = go.Figure(go.Sunburst(labels=data[LABEL], parents=data[PARENT], values=values,
                                    ids=data[IDS],
                                    hoverinfo='label+text', maxdepth=max_depth,
                                    branchvalues=branch_values,
                                    hovertext=get_hover_fig_text(data, TOPOLOGY_A, ref_base),
                                    marker=dict(colors=data[COUNT], colorscale=colorscale,
                                                cmin=c_min, cmax=c_max, cmid=c_mid, showscale=True,
                                                colorbar=dict(title=dict(text=colorbar_legend)))))
        fig.update_layout(title=dict(text=title, x=0.5, xanchor='center'))

    elif analysis == ENRICHMENT_A:
        data, significant = get_data_enrichment_analysis(data, ref_classes_abundance, test)
        data = data_cut_root(data, root_cut)
        fig = make_subplots(rows=1, cols=2,
                            column_widths=[0.3, 0.7],
                            vertical_spacing=0.03,
                            subplot_titles=(table_title, title),
                            specs=[[{'type': 'table'}, {'type': 'sunburst'}]])

        fig.add_trace(go.Sunburst(labels=data[LABEL], parents=data[PARENT],
                                  values=values, ids=data[IDS],
                                  hovertext=get_hover_fig_text(data, ENRICHMENT_A, ref_base),
                                  hoverinfo='label+text', maxdepth=max_depth,
                                  branchvalues=branch_values,
                                  marker=dict(colors=data[PVAL], colorscale=colorscale,
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
    if output is not None:
        fig.write_html(f'{output}.html')
    return fig


def get_hover_fig_text(data, analysis, ref_base):
    if analysis == ENRICHMENT_A:
        return [f'P value: {10 ** (-data[PVAL][i])}<br>'
                f'{COUNT}: <b>{data[COUNT][i]}</b><br>'
                f'{REF_COUNT}: {data[REF_COUNT][i]}<br>'
                f'{PROP}: <b>{round(data[PROP][i] * 100, 2)}%</b><br>'
                f'{REF_PROP}: {round(data[REF_PROP][i] * 100, 2)}%<br>'
                f'{IDS}: {data[ONTO_ID][i]}'
                if data[PVAL][i] > 0 else
                f'P value: {10 ** data[PVAL][i]}<br>'
                f'{COUNT}: <b>{data[COUNT][i]}</b><br>'
                f'{REF_COUNT}: {data[REF_COUNT][i]}<br>'
                f'{PROP}: <b>{round(data[PROP][i] * 100, 2)}%</b><br>'
                f'{REF_PROP}: {round(data[REF_PROP][i] * 100, 2)}%<br>'
                f'{IDS}: {data[ONTO_ID][i]}'
                for i in range(len(data[PVAL]))]
    elif analysis == TOPOLOGY_A:
        if ref_base:
            return [f'{COUNT}: <b>{data[COUNT][i]}</b><br>'
                    f'{REF_COUNT}: {data[REF_COUNT][i]}<br>'
                    f'{PROP}: <b>{round(data[PROP][i] * 100, 2)}%</b><br>'
                    f'{REF_PROP}: {round(data[REF_PROP][i] * 100, 2)}%<br>'
                    f'{IDS}: {data[ONTO_ID][i]}'
                    for i in range(len(data[PROP]))]
        else:
            return [f'{COUNT}: <b>{data[COUNT][i]}</b><br>'
                    f'{PROP}: <b>{round(data[PROP][i] * 100, 2)}%</b><br>'
                    f'{IDS}: {data[ONTO_ID][i]}'
                    for i in range(len(data[PROP]))]
