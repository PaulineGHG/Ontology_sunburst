from typing import List, Dict
from dash import Dash, dcc, html, Input, Output
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import scipy.stats as stats

# CONSTANTS
# ==================================================================================================
BINOMIAL_TEST = 'Binomial'
HYPERGEO_TEST = 'Hypergeometric'

# Keys
# ----
IDS = 'ID'
PARENT = 'Parent'
LABEL = 'Label'
COUNT = 'Count'
PROP = 'Proportion'
R_PROP = 'Relative_proportion'
PROP_DIF = 'Proportion_difference'
PVAL = 'Pvalue'


# FUNCTIONS
# ==================================================================================================
def get_fig_parameters(classes_abondance: Dict[str, int], parent_dict: Dict[str, List[str]],
                       children_dict: Dict[str, List[str]], root_item, full: bool = True,
                       names: Dict = None) -> Dict[str, List]:
    """ Returns a dictionary of parameters to create the sunburst figure.

    Parameters
    ----------
    classes_abondance: Dict[str, int]
        Dictionary associating for each class the number of metabolites found belonging to the class.
    parent_dict: Dict[str, List[str]]
        Dictionary associating for each class, its parents classes
    children_dict: Dict[str, List[str]]
        Dictionary associating for each class, its children classes
    full: bool (default=True)
        True for a full figure with class duplication if a class has +1 parents.
        False to have a reduced vew with all label appearing only once (1 random parent chosen)
    root_item: str
    names: dict

    Returns
    -------
    Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
    """
    data = {IDS: list(),
            PARENT: list(),
            LABEL: list(),
            COUNT: list()}
    for c_label, c_abundance in classes_abondance.items():
        if c_label != root_item:
            c_parents = parent_dict[c_label]
            if names is None:
                data = add_value_data(data, c_label, c_label, c_abundance, c_parents[0])
            else:
                data = add_value_data(data, c_label, names[c_label], c_abundance, c_parents[0])
            if len(c_parents) > 1 and full:
                for p in c_parents[1:]:
                    suffix = '__' + p
                    c_id = c_label + suffix
                    data = add_value_data(data, c_id, c_label, c_abundance, p)
                    if c_label in children_dict.keys():
                        c_children = children_dict[c_label]
                        for c in c_children:
                            data = add_children(data, suffix, c, c_id, classes_abondance,
                                                children_dict)
        else:
            data = add_value_data(data, c_label, c_label, c_abundance, '')

    return data


def add_value_data(data: Dict[str, List], m_id: str, label: str, value: int, parent: str) -> \
        Dict[str, List]:
    """ Fill the data dictionary for a metabolite class.

    Parameters
    ----------
    data: Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
    m_id: str
        ID of the metabolite class to add
    label: str
        Label (name) of the metabolite class to add
    value: int
        Abundance value of the metabolite class to add
    parent: str
        Parent metabolite class of the metabolite class to add

    Returns
    -------
    Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
    """
    data[IDS].append(m_id)
    data[LABEL].append(label)
    data[PARENT].append(parent)
    data[COUNT].append(value)
    return data


def add_children(data: Dict[str, List], origin: str, child: str, parent: str,
                 classes_abondance: Dict[str, int], children_dict: Dict[str, List[str]]) \
        -> Dict[str, List]:
    """ Add recursively all children of a given class to the data dictionary.

    Parameters
    ----------
    data: Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
    origin: str
        Origin of propagation : parent class of parent
    child: str
        Child metabolite class
    parent: str
        Parent metabolite class
    classes_abondance: Dict[str, int]
        Dictionary associating for each class the number of metabolites found belonging to the class.
    children_dict: Dict[str, List[str]]
        Dictionary associating for each class, its children classes

    Returns
    -------
    Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
    """
    if child in classes_abondance.keys():
        data[IDS].append(child + origin)
        data[LABEL].append(child)
        data[PARENT].append(parent)
        data[COUNT].append(classes_abondance[child])
        if child in children_dict.keys():
            origin_2 = origin + '__' + child
            cs = children_dict[child]
            for c in cs:
                add_children(data, origin_2, c, child + origin, classes_abondance, children_dict)
    return data


def get_data_proportion(data: Dict[str, List], total: bool) -> Dict[str, List]:
    """ Add a proportion value for color.

    Parameters
    ----------
    data: Dict[str, List]
    total: bool

    Returns
    -------
    Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
            - colors value : Proportion (0 < float <= 1)
            - branch proportion : Relative_prop
    """
    # Get total proportion
    max_abondance = np.max(data[COUNT])
    data[PROP] = [x / max_abondance for x in data[COUNT]]
    # Get proportion relative to +1 parent proportion for total branch value
    if total:
        data[R_PROP] = [x for x in data[PROP]]
        p = ''
        data = get_relative_prop(data, p)
    return data


def get_relative_prop(data, p):
    if p == '':
        prop_p = 1.0
        count_p = max(data[COUNT])
    else:
        prop_p = data[R_PROP][data[IDS].index(p)]
        count_p = data[COUNT][data[IDS].index(p)]
    index_p = [i for i, v in enumerate(data[PARENT]) if v == p]
    c_p = [data[IDS][i] for i in index_p]
    count_c_p = [data[COUNT][i] for i in index_p]
    if sum(count_c_p) > count_p:
        total = sum(count_c_p)
    else:
        total = count_p
    tot_prop = []
    for i, c in enumerate(c_p):
        prop = round((count_c_p[i] / total) * prop_p, 15)
        data[R_PROP][data[IDS].index(c)] = prop
        tot_prop.append(prop)
    # Correct prop to not exceed 100%
    while sum(tot_prop) > prop_p:
        tot_prop[-1] -= 0.000000000000001
    data[R_PROP][data[IDS].index(c)] = prop
    for c in c_p:
        if c in data[PARENT]:
            data = get_relative_prop(data, c)
    return data


def get_data_prop_diff(data, b_classes_abundance):
    i_max_abondance = np.max(data[COUNT])
    i_prop = [x / i_max_abondance for x in data[COUNT]]
    b_max_abondance = np.max(list(b_classes_abundance.values()))
    b_prop = [b_classes_abundance[x] / b_max_abondance if x in b_classes_abundance.keys() else 0 for
              x in data[LABEL]]
    diff = [i_prop[i] - b_prop[i] for i in range(len(i_prop))]
    data[PROP_DIF] = diff
    return data


def get_data_enrichment_analysis(data, b_classes_abundance, test, names):
    M = np.max(list(b_classes_abundance.values()))
    if names:
        m_list = [b_classes_abundance[x] if x in b_classes_abundance.keys() else 0 for x in
                  data[IDS]]
    else:
        m_list = [b_classes_abundance[x] if x in b_classes_abundance.keys() else 0 for x in
                  data[LABEL]]
    N = np.max(data[COUNT])
    n_list = data[COUNT]
    data[PVAL] = list()
    nb_classes = len(set(data[IDS]))
    significant_representation = dict()
    for i in range(len(m_list)):
        # Binomial Test
        if test == BINOMIAL_TEST:
            p_val = stats.binomtest(n_list[i], N, m_list[i] / M, alternative='two-sided').pvalue
        # Hypergeometric Test
        elif test == HYPERGEO_TEST:
            p_val = stats.hypergeom.sf(n_list[i] - 1, M, m_list[i], N)
        else:
            raise ValueError(f'test parameter must be in : {[BINOMIAL_TEST, HYPERGEO_TEST]}')
        if ((n_list[i] / N) - (m_list[i] / M)) > 0:
            data[PVAL].append(-np.log10(p_val))
        else:
            data[PVAL].append(np.log10(p_val))
        if p_val < 0.05 / nb_classes:
            significant_representation[data[IDS][i]] = p_val.round(10)
    significant_representation = dict(
        sorted(significant_representation.items(), key=lambda item: item[1]))
    return data, significant_representation


def generate_sunburst_fig(data: Dict[str, List[str or int or float]], output: str = None,
                          sb_type: str = 'proportion', b_classes_abond=None, test=BINOMIAL_TEST,
                          names: bool = False, total: bool = True):
    """ Generate a Sunburst figure and save it to output path.

    Parameters
    ----------
    data: Dict[str, List[str or int or float]]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
            - colors value : Proportion (0 < float <= 1)
    output: str
        Path to output to save the figure without extension
    sb_type: str
    b_classes_abond
    test
    names
    total
    """
    if total:
        branch_values = 'total'
        values = data[R_PROP]
    else:
        branch_values = 'remainder'
        values = data[COUNT]
    if sb_type == 'proportion':
        fig = go.Figure(go.Sunburst(labels=data[LABEL], parents=data[PARENT], values=values,
                                    ids=data[IDS],
                                    hoverinfo='label+text', maxdepth=7,
                                    branchvalues=branch_values,
                                    hovertext=[f'Count: {data[COUNT][i]}<br>'
                                               f'Proportion: {round(data[PROP][i]*100, 2)}%<br>'
                                               f'ID: {data[IDS][i]}'
                                               for i in range(len(data[PROP]))],
                                    marker=dict(colors=data[COUNT],
                                                colorscale=px.colors.diverging.curl_r,
                                                cmid=0.5 * max(data[COUNT]), showscale=True)))
    elif sb_type == 'comparison':
        data, signif = get_data_enrichment_analysis(data, b_classes_abond, test, names)
        m = np.mean(data[PROP_DIF])
        fig = make_subplots(rows=1, cols=2,
                            column_widths=[0.3, 0.7],
                            vertical_spacing=0.03,
                            subplot_titles=('Significant p-values',
                                            'Metabolites classes enrichment representation'),
                            specs=[[{'type': 'table'}, {'type': 'sunburst'}]])

        fig.add_trace(go.Sunburst(labels=data[LABEL], parents=data[PARENT],
                                  values=values, ids=data[IDS],
                                  hovertext=[f'P value: {10 ** (-data[PVAL][i])}<br>'
                                             f'Count: {data[COUNT][i]}<br>'
                                             f'Proportion: {round(data[PROP][i]*100, 2)}%<br>'
                                             f'ID: {data[IDS][i]}'
                                             if data[PVAL][i] > 0 else
                                             f'P value: {10 ** data[PVAL][i]}<br>'
                                             f'Count: {data[COUNT][i]}<br>'
                                             f'Proportion: {round(data[PROP][i]*100, 2)}%<br>'
                                             f'ID: {data[IDS][i]}'
                                             for i in range(len(data[PVAL]))],
                                  hoverinfo='label+text', maxdepth=7,
                                  branchvalues=branch_values,
                                  marker=dict(colors=data[PVAL],
                                              colorscale=px.colors.diverging.RdBu,
                                              cmid=0, cmax=10.0, cmin=-10.0, showscale=True)),
                      row=1, col=2)

        fig.add_trace(go.Table(header=dict(values=['Metabolite', f'{test} test P-value'],
                                           fill=dict(color='#666666'), height=40,
                                           font=dict(size=20)),
                               cells=dict(values=[list(signif.keys()), list(signif.values())],
                                          fill=dict(color='#777777'), height=35,
                                          font=dict(size=16))),
                      row=1, col=1)
    else:
        raise ValueError('Wrong type input')
    fig.update_layout(paper_bgcolor="#888888", font_color='#111111', font_size=20)
    fig.update_annotations(font_size=28)
    if output is not None:
        fig.write_html(f'{output}.html')
    return fig

# def generate_dash_interactive_fig(data, output):
#     app = Dash(__name__)
#     values = ['Proportion', 'Proportion Difference']
#     color_scale = {'Proportion': px.colors.sequential.Mint,
#                    'Proportion Difference': px.colors.diverging.Armyrose_r}
#     midpoint = {'Proportion': 0.5,
#                 'Proportion Difference': np.mean(data['Proportion Difference'])}
#     app.layout = html.Div([html.P("Select value:"),
#                            dcc.Dropdown(id='dropdown',
#                                         options=values,
#                                         value='Proportion'),
#                            dcc.Graph(id="graph")])
#
#     @app.callback(
#         Output("graph", "figure"),
#         Input("dropdown", "value"))
#     def change_measure(value):
#         fig = px.sunburst(data, names='Label', parents='Parent', values='Count', ids='ID', color=value,
#                           template='presentation', color_continuous_scale=color_scale[value], height=900,
#                           color_continuous_midpoint=midpoint[value])
#         fig.update_layout(paper_bgcolor="#373f48", font_color='#ffffff')
#         return fig
#     app.run_server(debug=True)
