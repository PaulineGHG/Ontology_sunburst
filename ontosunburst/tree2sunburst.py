import os.path
import difflib
import inspect
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from ontosunburst.dag2tree import *
from ontosunburst.onto2dag import *

# ==================================================================================================
# CONSTANTS
# ==================================================================================================

# COLORIZATION
Colorization: TypeAlias = Literal['w', 'p', 'e', 'd']
COL_W, COL_P, COL_E, COL_D = get_args(Colorization)
COLOR_DATA = {COL_W: CW, COL_P: PROP, COL_E: E_L10PV, COL_D: DIFF}
COLOR_DATA_REF = {COL_W: R_CW, COL_P: R_PROP, COL_E: E_L10PV, COL_D: DIFF}

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
KWARGS = [C_MIN, C_MAX, C_MID, MAX_DEPTH, COLORSCALE, TITLE, COLORBAR_LEGEND, BG_COLOR, FONT_COLOR,
          FONT_SIZE, TABLE_TITLE, TABLE_LEGEND, TABLE_COLOR]
KWARGS_TYPE = {C_MIN: float, C_MAX: float, C_MID: float, MAX_DEPTH: int, COLORSCALE: str,
               TITLE: str, COLORBAR_LEGEND: str, BG_COLOR: str, FONT_COLOR: str, FONT_SIZE: int,
               TABLE_TITLE: str, TABLE_LEGEND: str, TABLE_COLOR: str}


# ==================================================================================================
# FUNCTIONS
# ==================================================================================================
def generate_sunburst_fig(tree: InducedTree, output: str, colorization: str, sector_sized: bool,
                          ref_set: bool, ref_base: bool, write_fig: bool = True, **kwargs) \
        -> go.Figure:
    c_min, c_max, c_mid, max_depth, colorscale, title, colorbar_legend, background_color, \
        font_color, font_size, table_title, table_legend, table_color = \
        get_fig_kwargs(output, colorization, **kwargs)
    data = get_tree_data(tree)
    fig = go.Figure(go.Sunburst(labels=data[LAB], parents=data[PARENT],
                                values=data[RELAT_P], ids=data[ID],
                                hoverinfo='label+text', maxdepth=max_depth,
                                branchvalues='total',
                                hovertext=get_hover_fig_text(data, ref_set),
                                marker=dict(colors=data[CW], colorscale=colorscale,
                                            cmin=c_min, cmax=c_max, cmid=c_mid, showscale=True,
                                            colorbar=dict(title=dict(text=colorbar_legend)))))

    fig.update_layout(title=dict(text=title, x=0.5, xanchor='center'))

    fig.update_layout(paper_bgcolor=background_color, font_color=font_color, font_size=font_size)
    fig.update_annotations(font_size=font_size * 1.5)
    if write_fig:
        fig.write_html(f'{output}.html')
        # write_tsv_output(data, f'{output}.tsv')
    return fig


def get_tree_data(tree: InducedTree):
    root_node = tree.nodes[(tree.sub_dag.root, 0)]
    data = {x: [] for x in root_node.get_arguments().keys()}

    for node in tree.nodes.values():
        node_args = node.get_arguments()
        for arg, val in node_args.items():
            data[arg].append(val)
    print(data)
    return data


# Figure creation
# --------------------------------------------------------------------------------------------------
def get_fig_kwargs(output: str, colorization: str, **kwargs):
    """ Generate a Sunburst figure and save it to output path.

        Parameters
        ----------
        output: str (optional, default=None)
            Path to output to save the figure without extension

        """
    check_kwargs(**kwargs)
    def_colorscale = {COL_W: 'Viridis',
                      COL_P: 'Viridis',
                      COL_E: 'RdBu',
                      COL_D: 'RdBu'}
    def_titles = {COL_W: f'{os.path.basename(output)} : Cumulative weight',
                  COL_P: f'{os.path.basename(output)} : Proportion',
                  COL_E: f'{os.path.basename(output)} : Enrichment representation',
                  COL_D: f'{os.path.basename(output)} : Proportion difference'}
    def_colorbar = {COL_W: 'Weight',
                    COL_P: '%age',
                    COL_E: 'Log10(p-value)',
                    COL_D: '%age diff'}
    def_c_min = {COL_W: 1, COL_P: 0, COL_E: -10, COL_D: -1}
    def_c_max = {COL_W: None, COL_P: 1,  COL_E: 10, COL_D: 1}
    def_c_mid = {COL_W: None, COL_P: 0.5, COL_E: 0, COL_D: 0}

    c_min = kwargs.get(C_MIN, def_c_min[colorization])
    c_max = kwargs.get(C_MAX, def_c_max[colorization])
    c_mid = kwargs.get(C_MID, def_c_mid[colorization])
    max_depth = kwargs.get(MAX_DEPTH, 7)
    colorscale = px.colors.get_colorscale(kwargs.get(COLORSCALE, def_colorscale[colorization]))
    title = kwargs.get(TITLE, def_titles[colorization])
    colorbar_legend = kwargs.get(COLORBAR_LEGEND, def_colorbar[colorization])
    background_color = kwargs.get(BG_COLOR, 'rgba(255, 255, 255, 0)')
    font_color = kwargs.get(FONT_COLOR, '#111111')
    font_size = kwargs.get(FONT_SIZE, 20)
    table_title = kwargs.get(TABLE_TITLE, 'Significant p-values')
    table_legend = kwargs.get(TABLE_LEGEND, 'IDs')
    table_color = kwargs.get(TABLE_COLOR, '#666666')

    return c_min, c_max, c_mid, max_depth, colorscale, title, colorbar_legend, background_color, \
        font_color, font_size, table_title, table_legend, table_color


def check_kwargs(**kwargs):
    close_matches = {x: difflib.get_close_matches(x, KWARGS, n=1, cutoff=0.5)[0] for x in kwargs
                     if difflib.get_close_matches(x, KWARGS, n=1, cutoff=0.5) and x not in KWARGS}
    for k in kwargs:
        if k not in KWARGS:
            if k in close_matches:
                print(f'Unknown kwarg "{k}", did you mean "{close_matches[k]}" ?')
            else:
                print(f'Unknown kwarg "{k}"')
        elif type(k) != KWARGS_TYPE[k]:
            print(f'"{k}" must be of type "{KWARGS_TYPE[k]}" not "{type(k)}"')


# def generate_sunburst_fig(data: TreeData, output: str, analysis: str = TOPOLOGY_A,
#                           test=BINOMIAL_TEST, significant: Dict[str, float] = None,
#                           ref_set: bool = True, write_fig: bool = True, **kwargs) -> go.Figure:
#     """ Generate a Sunburst figure and save it to output path.
#
#     Parameters
#     ----------
#     data: TreeData
#         DataTable of figure parameters
#         (sectors id, label, parent, count, proportion, p-value, ...)
#     output: str
#         Path to output to save the figure without extension
#     analysis: str (optional, default=topology)
#         Analysis mode : topology or enrichment
#     test: str (optional, default=Binomial)
#         Type of test for enrichment analysis : binomial or hypergeometric
#     significant: Dict[str, float]
#         Dictionary of significant p-value {ontology-id: p-value}
#     ref_set: bool (optional, default=True)
#         True if a reference set is present, False otherwise. If true will show reference set values
#         in the hover text of sectors.
#     write_fig: bool (optional, default=True)
#         True to write the html figure, False to only return figure
#     **kwargs
#         Keyword args: c_min, c_max, c_mid, max_depth, colorscale, title, colorbar_legend, bg_color,
#         font_color, font_size, table_title, table_legend, table_color
#
#     Returns
#     -------
#     go.Figure
#         Sunburst figure generated.
#     """
#     c_min, c_max, c_mid, max_depth, colorscale, title, colorbar_legend, background_color, \
#         font_color, font_size, table_title, table_legend, table_color = \
#         get_fig_kwargs(output, analysis, **kwargs)
#
#     if analysis == TOPOLOGY_A:
#         fig = go.Figure(go.Sunburst(labels=data.labels, parents=data.parents,
#                                     values=data.relative_prop, ids=data.ids,
#                                     hoverinfo='label+text', maxdepth=max_depth,
#                                     branchvalues='total',
#                                     hovertext=get_hover_fig_text(data, TOPOLOGY_A, ref_set),
#                                     marker=dict(colors=data.count, colorscale=colorscale,
#                                                 cmin=c_min, cmax=c_max, cmid=c_mid, showscale=True,
#                                                 colorbar=dict(title=dict(text=colorbar_legend)))))
#         fig.update_layout(title=dict(text=title, x=0.5, xanchor='center'))
#
#     elif analysis == ENRICHMENT_A:
#         fig = make_subplots(rows=1, cols=2,
#                             column_widths=[0.3, 0.7],
#                             vertical_spacing=0.03,
#                             subplot_titles=(table_title, title),
#                             specs=[[{'type': 'table'}, {'type': 'sunburst'}]])
#
#         fig.add_trace(go.Sunburst(labels=data.labels, parents=data.parents,
#                                   values=data.relative_prop, ids=data.ids,
#                                   hovertext=get_hover_fig_text(data, ENRICHMENT_A, ref_set),
#                                   hoverinfo='label+text', maxdepth=max_depth,
#                                   branchvalues='total',
#                                   marker=dict(colors=data.p_val, colorscale=colorscale,
#                                               cmid=c_mid, cmax=c_max, cmin=c_min, showscale=True,
#                                               colorbar=dict(title=dict(text=colorbar_legend)))),
#                       row=1, col=2)
#
#         fig.add_trace(go.Table(header=dict(values=[table_legend, f'{test} test P-value'],
#                                            fill=dict(color=table_color), height=40,
#                                            font=dict(size=font_size)),
#                                cells=dict(values=[list(significant.keys()),
#                                                   list(significant.values())],
#                                           fill=dict(color=table_color), height=35,
#                                           font=dict(size=font_size * 0.80))),
#                       row=1, col=1)
#     else:
#         raise ValueError('Wrong type input')
#     fig.update_layout(paper_bgcolor=background_color, font_color=font_color, font_size=font_size)
#     fig.update_annotations(font_size=font_size * 1.5)
#     if write_fig:
#         fig.write_html(f'{output}.html')
#         write_tsv_output(data, f'{output}.tsv')
#     return fig


def get_hover_fig_text(data: Dict, ref_set: bool) \
        -> List[str]:
    """

    Parameters
    ----------
    data
    ref_set

    Returns
    -------

    """
    data_len = len(data[CW])
    if ref_set:
        return [f'{CW}: <b>{data[CW][i]}</b><br>'
                f'{R_CW}: {data[R_CW][i]}<br>'
                f'{PROP}: <b>{round(data[PROP][i] * 100, 2)}%</b><br>'
                f'{R_PROP}: {round(data[R_PROP][i] * 100, 2)}%<br>'
                f'{EW}: {data[EW][i]}<br>'
                f'{R_EW}: {data[R_EW][i]}<br>'
                f'{E_PV}: {data[E_PV][i]}<br>'
                f'{DIFF}: {data[DIFF][i]}<br>'
                f'{O_ID}: {data[O_ID][i]}'
                for i in range(data_len)]
    else:
        return [f'{CW}: <b>{data[CW][i]}</b><br>'
                f'{PROP}: <b>{round(data[PROP][i] * 100, 2)}%</b><br>'
                f'{EW}: <b>{data[EW][i]}</b><br>'
                f'{O_ID}: {data[O_ID][i]}'
                for i in range(data_len)]
#
#
# def write_tsv_output(data, output):
#     d_data = data.get_data_dict()
#     d_data_id = dict()
#     for i in range(data.len):
#         if d_data[ONTO_ID][i] not in d_data_id:
#             d_data_id[d_data[ONTO_ID][i]] = {'Parents ids': [],
#                                              'Parents labels': [],
#                                              LABEL: d_data[LABEL][i].replace(' ...', ''),
#                                              WEIGHT: d_data[WEIGHT][i],
#                                              REF_WEIGHT: d_data[REF_WEIGHT][i],
#                                              PROP: d_data[PROP][i],
#                                              REF_PROP: d_data[REF_PROP][i],
#                                              PVAL: d_data[PVAL][i]}
#         parent_id = d_data[PARENT][i]
#         if parent_id in d_data[IDS]:
#             p_index = d_data[IDS].index(parent_id)
#             parent_onto_id = d_data[ONTO_ID][p_index]
#             parent_label = d_data[LABEL][p_index].replace(' ...', '')
#             d_data_id[d_data[ONTO_ID][i]]['Parents ids'].append(parent_onto_id)
#             d_data_id[d_data[ONTO_ID][i]]['Parents labels'].append(parent_label)
#         else:
#             d_data_id[d_data[ONTO_ID][i]]['Parents ids'].append(parent_id)
#             d_data_id[d_data[ONTO_ID][i]]['Parents labels'].append(parent_id)
#
#     print(len(d_data_id))
#
#     with open(output, 'w') as f:
#         f.write(f'{ONTO_ID}\t{LABEL}\tParents ids\tParents labels\t'
#                 f'{WEIGHT}\t{REF_WEIGHT}\t{PROP}\t{REF_PROP}\t{PVAL}\n')
#         for onto_id, val in d_data_id.items():
#             f.write(f'{onto_id}\t{val[LABEL]}\t{" ; ".join(val["Parents ids"])}\t'
#                     f'{" ; ".join(val["Parents labels"])}\t{val[WEIGHT]}\t{val[REF_WEIGHT]}\t'
#                     f'{val[PROP]}\t{val[REF_PROP]}\t{val[PVAL]}\n')

