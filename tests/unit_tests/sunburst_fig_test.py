import unittest
from unittest.mock import patch
import io
import sys
from functools import wraps

import plotly.graph_objs
from numpy import nan
from ontosunburst.sunburst_fig import *
from ontosunburst.ontology import *

"""
Tests manually good file creation.
No automatic tests integrated.
"""

# ==================================================================================================
# GLOBAL
# ==================================================================================================

# GENERAL DICT ONTO (METACYC, KEGG)
# --------------------------------------------------------------------------------------------------
ENRICH_DATA = {IDS: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
               ONTO_ID: ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09'],
               PARENT: ['', 0, 0, 0, 0, 1, 1, 1, 2, 2],
               LABEL: ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'],
               COUNT: [50, 5, 25, 20, 1, 5, nan, nan, 1, 1],
               REF_COUNT: [100, 40, 30, 20, 10, 20, 5, 1, 1, 3],
               PROP: [1, 0.1, 0.5, 0.4, 0.02, 0.1, nan, nan, 0.02, 0.02],
               REF_PROP: [1, 0.4, 0.3, 0.2, 0.1, 0.2, 0.05, 0.01, 0.01, 0.03],
               RELAT_PROP: [1000000, 400000.0, 300000.0, 200000.0, 100000.0, 200000.0, 50000.0,
                            10000.0, 10000.0, 30000.0]}
p_val = [0.0, -5.420266413988895, 2.509702991379166, 2.948803113091024,
         -1.2341542222355069, -1.103304935668835, nan, nan, 0.4034095751193356, 0.0]
ENRICH_REF_AB = {'00': 100, '01': 40, '02': 30, '03': 20, '04': 10, '05': 20, '06': 5, '07': 1,
                 '08': 1, '09': 3}


# ==================================================================================================
# FUNCTIONS UTILS
# ==================================================================================================
def dicts_with_sorted_lists_equal(dict1, dict2):
    if dict1.keys() != dict2.keys():
        return False
    for key in dict1:
        if sorted(dict1[key]) != sorted(dict2[key]):
            return False
    return True


def test_for(func):
    def decorator(test_func):
        @wraps(test_func)
        def wrapper(*args, **kwargs):
            return test_func(*args, **kwargs)

        wrapper._test_for = func
        return wrapper

    return decorator


class DualWriter(io.StringIO):
    def __init__(self, original_stdout):
        super().__init__()
        self.original_stdout = original_stdout

    def write(self, s):
        super().write(s)
        self.original_stdout.write(s)


# ==================================================================================================
# UNIT TESTS
# ==================================================================================================

# TEST
# --------------------------------------------------------------------------------------------------

class TestSunburstFigure(unittest.TestCase):
    @test_for(get_fig_kwargs)
    def test_get_kwargs_topology(self):
        c_min, c_max, c_mid, max_depth, colorscale, title, colorbar_legend, background_color, \
            font_color, font_size, table_title, table_legend, table_color = \
            get_fig_kwargs('out', TOPOLOGY_A)
        self.assertEqual(c_min, 1)
        self.assertEqual(c_max, None)
        self.assertEqual(c_mid, None)
        self.assertEqual(max_depth, 7)
        self.assertEqual(colorscale, px.colors.get_colorscale('Viridis'))
        self.assertEqual(title, 'out : Proportion of classes')
        self.assertEqual(colorbar_legend, 'Count')
        self.assertEqual(background_color, 'rgba(255, 255, 255, 0)')
        self.assertEqual(font_color, '#111111')
        self.assertEqual(font_size, 20)
        self.assertEqual(table_title, 'Significant p-values')
        self.assertEqual(table_legend, 'IDs')
        self.assertEqual(table_color, '#666666')

    @test_for(get_fig_kwargs)
    def test_get_kwargs_enrichment(self):
        c_min, c_max, c_mid, max_depth, colorscale, title, colorbar_legend, background_color, \
            font_color, font_size, table_title, table_legend, table_color = \
            get_fig_kwargs('out', ENRICHMENT_A)
        self.assertEqual(c_min, -10)
        self.assertEqual(c_max, 10)
        self.assertEqual(c_mid, 0)
        self.assertEqual(max_depth, 7)
        self.assertEqual(colorscale, px.colors.get_colorscale('RdBu'))
        self.assertEqual(title, 'out : Classes enrichment representation')
        self.assertEqual(colorbar_legend, 'Log10(p-value)')
        self.assertEqual(background_color, 'rgba(255, 255, 255, 0)')
        self.assertEqual(font_color, '#111111')
        self.assertEqual(font_size, 20)
        self.assertEqual(table_title, 'Significant p-values')
        self.assertEqual(table_legend, 'IDs')
        self.assertEqual(table_color, '#666666')

    @test_for(get_fig_kwargs)
    def test_get_kwargs(self):
        c_min, c_max, c_mid, max_depth, colorscale, title, colorbar_legend, background_color, \
            font_color, font_size, table_title, table_legend, table_color = \
            get_fig_kwargs('out', TOPOLOGY_A, c_min=4, c_max=8, c_mid=6, max_depth=10,
                           colorscale='Twilight', title='My title', colorbar_legend='Total',
                           bg_color='#000000', font_color='#ffffff', font_size=24,
                           table_title='Significant classes', table_legend='Name',
                           table_color='#222222')
        self.assertEqual(c_min, 4)
        self.assertEqual(c_max, 8)
        self.assertEqual(c_mid, 6)
        self.assertEqual(max_depth, 10)
        self.assertEqual(colorscale, px.colors.get_colorscale('Twilight'))
        self.assertEqual(title, 'My title')
        self.assertEqual(colorbar_legend, 'Total')
        self.assertEqual(background_color, '#000000')
        self.assertEqual(font_color, '#ffffff')
        self.assertEqual(font_size, 24)
        self.assertEqual(table_title, 'Significant classes')
        self.assertEqual(table_legend, 'Name')
        self.assertEqual(table_color, '#222222')

    @test_for(get_hover_fig_text)
    def test_get_hover_fig_text_enrich_ref(self):
        text_list = get_hover_fig_text(ENRICH_DATA, ENRICHMENT_A, True)
        self.assertEqual(len(text_list), 10)
        self.assertEqual(text_list[0], 'P value: 1.0<br>Count: <b>50</b>'
                                       '<br>Reference count: 100<br>Proportion: <b>100%</b>'
                                       '<br>Reference proportion: 100%<br>ID: 00')
        self.assertEqual(text_list[5], 'P value: 0.07883064215278136<br>Count: <b>5</b>'
                                       '<br>Reference count: 20<br>Proportion: <b>10.0%</b>'
                                       '<br>Reference proportion: 20.0%<br>ID: 05')

    @test_for(get_hover_fig_text)
    def test_get_hover_fig_text_enrich_no_ref(self):
        text_list = get_hover_fig_text(ENRICH_DATA, ENRICHMENT_A, False)
        self.assertEqual(len(text_list), 10)
        self.assertEqual(text_list[0], 'P value: 1.0<br>Count: <b>50</b>'
                                       '<br>Reference count: 100<br>Proportion: <b>100%</b>'
                                       '<br>Reference proportion: 100%<br>ID: 00')
        self.assertEqual(text_list[5], 'P value: 0.07883064215278136<br>Count: <b>5</b>'
                                       '<br>Reference count: 20<br>Proportion: <b>10.0%</b>'
                                       '<br>Reference proportion: 20.0%<br>ID: 05')

    @test_for(get_hover_fig_text)
    def test_get_hover_fig_text_topology_ref(self):
        text_list = get_hover_fig_text(ENRICH_DATA, TOPOLOGY_A, True)
        self.assertEqual(len(text_list), 10)
        self.assertEqual(text_list[0], 'Count: <b>50</b><br>Reference count: 100'
                                       '<br>Proportion: <b>100%</b>'
                                       '<br>Reference proportion: 100%<br>ID: 00')
        self.assertEqual(text_list[5], 'Count: <b>5</b><br>Reference count: 20'
                                       '<br>Proportion: <b>10.0%</b>'
                                       '<br>Reference proportion: 20.0%<br>ID: 05')

    @test_for(get_hover_fig_text)
    def test_get_hover_fig_text_topology_no_ref(self):
        text_list = get_hover_fig_text(ENRICH_DATA, TOPOLOGY_A, False)
        self.assertEqual(len(text_list), 10)
        self.assertEqual(text_list[0], 'Count: <b>50</b><br>Proportion: <b>100%</b><br>ID: 00')
        self.assertEqual(text_list[5], 'Count: <b>5</b><br>Proportion: <b>10.0%</b><br>ID: 05')

    @test_for(generate_sunburst_fig)
    def test_generate_sunburst_fig(self):
        fig = generate_sunburst_fig(ENRICH_DATA, 'out', analysis=ENRICHMENT_A,
                                    ref_classes_abundance=ENRICH_REF_AB, test=HYPERGEO_TEST)
        fig = fig.to_dict()
        w_fig = {'data': [{'branchvalues': 'total', 'hoverinfo': 'label+text',
                           'hovertext': ['P value: 7.263166523971595e-10<br>Count: <b>5</b>'
                                         '<br>Reference count: 40<br>Proportion: <b>10.0%</b>'
                                         '<br>Reference proportion: 40.0%<br>ID: 01',
                                         'P value: 2.029502412840083e-05<br>Count: <b>25</b>'
                                         '<br>Reference count: 30<br>Proportion: <b>50.0%</b>'
                                         '<br>Reference proportion: 30.0%<br>ID: 02',
                                         'P value: 1.7586072571039989e-07<br>Count: <b>20</b>'
                                         '<br>Reference count: 20<br>Proportion: <b>40.0%</b>'
                                         '<br>Reference proportion: 20.0%<br>ID: 03',
                                         'P value: 0.015660489896045526<br>Count: <b>1</b>'
                                         '<br>Reference count: 10<br>Proportion: <b>2.0%</b>'
                                         '<br>Reference proportion: 10.0%<br>ID: 04',
                                         'P value: 0.022834981011415515<br>Count: <b>5</b>'
                                         '<br>Reference count: 20<br>Proportion: <b>10.0%</b>'
                                         '<br>Reference proportion: 20.0%<br>ID: 05',
                                         'P value: nan<br>Count: <b>nan</b>'
                                         '<br>Reference count: 5<br>Proportion: <b>nan%</b>'
                                         '<br>Reference proportion: 5.0%<br>ID: 06',
                                         'P value: nan<br>Count: <b>nan</b>'
                                         '<br>Reference count: 1<br>Proportion: <b>nan%</b>'
                                         '<br>Reference proportion: 1.0%<br>ID: 07',
                                         'P value: 1.0<br>Count: <b>1</b>'
                                         '<br>Reference count: 1<br>Proportion: <b>2.0%</b>'
                                         '<br>Reference proportion: 1.0%<br>ID: 08',
                                         'P value: 0.9999999999999997<br>Count: <b>1</b>'
                                         '<br>Reference count: 3<br>Proportion: <b>2.0%</b>'
                                         '<br>Reference proportion: 3.0%<br>ID: 09'],
                           'ids': [1, 2, 3, 4, 5, 6, 7, 8, 9],
                           'labels': ['1', '2', '3', '4', '5', '6', '7', '8', '9'],
                           'marker': {'cmax': 10, 'cmid': 0, 'cmin': -10,
                                      'colorbar': {'title': {'text': 'Log10(p-value)'}},
                                      'colors': [-9.138873998573988, 4.692610428021241,
                                                 6.754831139005899, -1.8051946563380086,
                                                 -1.6413993451973743, nan, nan, -0.0,
                                                 -1.4464911998299308e-16],
                                      'colorscale': [[0.0, 'rgb(103,0,31)'],
                                                     [0.1, 'rgb(178,24,43)'],
                                                     [0.2, 'rgb(214,96,77)'],
                                                     [0.30000000000000004, 'rgb(244,165,130)'],
                                                     [0.4, 'rgb(253,219,199)'],
                                                     [0.5, 'rgb(247,247,247)'],
                                                     [0.6000000000000001, 'rgb(209,229,240)'],
                                                     [0.7000000000000001, 'rgb(146,197,222)'],
                                                     [0.8, 'rgb(67,147,195)'],
                                                     [0.9, 'rgb(33,102,172)'],
                                                     [1.0, 'rgb(5,48,97)']],
                                      'showscale': True},
                           'maxdepth': 7, 'parents': [0, 0, 0, 0, 1, 1, 1, 2, 2],
                           'values': [400000.0, 300000.0, 200000.0, 100000.0, 200000.0, 50000.0,
                                      10000.0, 10000.0, 30000.0], 'type': 'sunburst',
                           'domain': {'x': [0.37, 1.0], 'y': [0.0, 1.0]}},
                          {'cells': {'fill': {'color': '#666666'}, 'font': {'size': 16.0},
                                     'height': 35, 'values': [['01', '03', '02'],
                                                              [7e-10, 1.759e-07, 2.0295e-05]]},
                           'header': {'fill': {'color': '#666666'}, 'font': {'size': 20},
                                      'height': 40,
                                      'values': ['IDs', 'hypergeometric test P-value']},
                           'type': 'table', 'domain': {'x': [0.0, 0.27], 'y': [0.0, 1.0]}}],
                 'layout': {'template': {'data': {'histogram2dcontour':
                                                      [{'type': 'histogram2dcontour',
                                                        'colorbar': {'outlinewidth': 0,
                                                                     'ticks': ''},
                                                        'colorscale': [[0.0, '#0d0887'],
                                                                       [0.1111111111111111,
                                                                        '#46039f'],
                                                                       [0.2222222222222222,
                                                                        '#7201a8'],
                                                                       [0.3333333333333333,
                                                                        '#9c179e'],
                                                                       [0.4444444444444444,
                                                                        '#bd3786'],
                                                                       [0.5555555555555556,
                                                                        '#d8576b'],
                                                                       [0.6666666666666666,
                                                                        '#ed7953'],
                                                                       [0.7777777777777778,
                                                                        '#fb9f3a'],
                                                                       [0.8888888888888888,
                                                                        '#fdca26'],
                                                                       [1.0, '#f0f921']]}],
                                                  'choropleth': [{'type': 'choropleth',
                                                                  'colorbar': {'outlinewidth': 0,
                                                                               'ticks': ''}}],
                                                  'histogram2d': [{'type': 'histogram2d',
                                                                   'colorbar': {'outlinewidth': 0,
                                                                                'ticks': ''},
                                                                   'colorscale': [[0.0, '#0d0887'],
                                                                                  [
                                                                                      0.1111111111111111,
                                                                                      '#46039f'], [
                                                                                      0.2222222222222222,
                                                                                      '#7201a8'], [
                                                                                      0.3333333333333333,
                                                                                      '#9c179e'], [
                                                                                      0.4444444444444444,
                                                                                      '#bd3786'], [
                                                                                      0.5555555555555556,
                                                                                      '#d8576b'], [
                                                                                      0.6666666666666666,
                                                                                      '#ed7953'], [
                                                                                      0.7777777777777778,
                                                                                      '#fb9f3a'], [
                                                                                      0.8888888888888888,
                                                                                      '#fdca26'],
                                                                                  [1.0,
                                                                                   '#f0f921']]}],
                                                  'heatmap': [{'type': 'heatmap',
                                                               'colorbar': {'outlinewidth': 0,
                                                                            'ticks': ''},
                                                               'colorscale': [[0.0, '#0d0887'],
                                                                              [0.1111111111111111,
                                                                               '#46039f'],
                                                                              [0.2222222222222222,
                                                                               '#7201a8'],
                                                                              [0.3333333333333333,
                                                                               '#9c179e'],
                                                                              [0.4444444444444444,
                                                                               '#bd3786'],
                                                                              [0.5555555555555556,
                                                                               '#d8576b'],
                                                                              [0.6666666666666666,
                                                                               '#ed7953'],
                                                                              [0.7777777777777778,
                                                                               '#fb9f3a'],
                                                                              [0.8888888888888888,
                                                                               '#fdca26'],
                                                                              [1.0, '#f0f921']]}],
                                                  'heatmapgl': [{'type': 'heatmapgl',
                                                                 'colorbar': {'outlinewidth': 0,
                                                                              'ticks': ''},
                                                                 'colorscale': [[0.0, '#0d0887'],
                                                                                [0.1111111111111111,
                                                                                 '#46039f'],
                                                                                [0.2222222222222222,
                                                                                 '#7201a8'],
                                                                                [0.3333333333333333,
                                                                                 '#9c179e'],
                                                                                [0.4444444444444444,
                                                                                 '#bd3786'],
                                                                                [0.5555555555555556,
                                                                                 '#d8576b'],
                                                                                [0.6666666666666666,
                                                                                 '#ed7953'],
                                                                                [0.7777777777777778,
                                                                                 '#fb9f3a'],
                                                                                [0.8888888888888888,
                                                                                 '#fdca26'],
                                                                                [1.0, '#f0f921']]}],
                                                  'contourcarpet': [{'type': 'contourcarpet',
                                                                     'colorbar': {'outlinewidth': 0,
                                                                                  'ticks': ''}}],
                                                  'contour': [{'type': 'contour',
                                                               'colorbar': {'outlinewidth': 0,
                                                                            'ticks': ''},
                                                               'colorscale': [[0.0, '#0d0887'],
                                                                              [0.1111111111111111,
                                                                               '#46039f'],
                                                                              [0.2222222222222222,
                                                                               '#7201a8'],
                                                                              [0.3333333333333333,
                                                                               '#9c179e'],
                                                                              [0.4444444444444444,
                                                                               '#bd3786'],
                                                                              [0.5555555555555556,
                                                                               '#d8576b'],
                                                                              [0.6666666666666666,
                                                                               '#ed7953'],
                                                                              [0.7777777777777778,
                                                                               '#fb9f3a'],
                                                                              [0.8888888888888888,
                                                                               '#fdca26'],
                                                                              [1.0, '#f0f921']]}],
                                                  'surface': [{'type': 'surface',
                                                               'colorbar': {'outlinewidth': 0,
                                                                            'ticks': ''},
                                                               'colorscale': [[0.0, '#0d0887'],
                                                                              [0.1111111111111111,
                                                                               '#46039f'],
                                                                              [0.2222222222222222,
                                                                               '#7201a8'],
                                                                              [0.3333333333333333,
                                                                               '#9c179e'],
                                                                              [0.4444444444444444,
                                                                               '#bd3786'],
                                                                              [0.5555555555555556,
                                                                               '#d8576b'],
                                                                              [0.6666666666666666,
                                                                               '#ed7953'],
                                                                              [0.7777777777777778,
                                                                               '#fb9f3a'],
                                                                              [0.8888888888888888,
                                                                               '#fdca26'],
                                                                              [1.0, '#f0f921']]}],
                                                  'mesh3d': [{'type': 'mesh3d',
                                                              'colorbar': {'outlinewidth': 0,
                                                                           'ticks': ''}}],
                                                  'scatter': [{'fillpattern': {
                                                      'fillmode': 'overlay', 'size': 10,
                                                      'solidity': 0.2}, 'type': 'scatter'}],
                                                  'parcoords': [{'type': 'parcoords', 'line': {
                                                      'colorbar': {'outlinewidth': 0,
                                                                   'ticks': ''}}}],
                                                  'scatterpolargl': [{'type': 'scatterpolargl',
                                                                      'marker': {'colorbar': {
                                                                          'outlinewidth': 0,
                                                                          'ticks': ''}}}], 'bar': [
                         {'error_x': {'color': '#2a3f5f'}, 'error_y': {'color': '#2a3f5f'},
                          'marker': {'line': {'color': '#E5ECF6', 'width': 0.5},
                                     'pattern': {'fillmode': 'overlay', 'size': 10,
                                                 'solidity': 0.2}}, 'type': 'bar'}], 'scattergeo': [
                         {'type': 'scattergeo',
                          'marker': {'colorbar': {'outlinewidth': 0, 'ticks': ''}}}],
                                                  'scatterpolar': [{'type': 'scatterpolar',
                                                                    'marker': {'colorbar': {
                                                                        'outlinewidth': 0,
                                                                        'ticks': ''}}}],
                                                  'histogram': [{'marker': {
                                                      'pattern': {'fillmode': 'overlay', 'size': 10,
                                                                  'solidity': 0.2}},
                                                                 'type': 'histogram'}],
                                                  'scattergl': [{'type': 'scattergl', 'marker': {
                                                      'colorbar': {'outlinewidth': 0,
                                                                   'ticks': ''}}}], 'scatter3d': [
                         {'type': 'scatter3d',
                          'line': {'colorbar': {'outlinewidth': 0, 'ticks': ''}},
                          'marker': {'colorbar': {'outlinewidth': 0, 'ticks': ''}}}],
                                                  'scattermapbox': [{'type': 'scattermapbox',
                                                                     'marker': {'colorbar': {
                                                                         'outlinewidth': 0,
                                                                         'ticks': ''}}}],
                                                  'scatterternary': [{'type': 'scatterternary',
                                                                      'marker': {'colorbar': {
                                                                          'outlinewidth': 0,
                                                                          'ticks': ''}}}],
                                                  'scattercarpet': [{'type': 'scattercarpet',
                                                                     'marker': {'colorbar': {
                                                                         'outlinewidth': 0,
                                                                         'ticks': ''}}}],
                                                  'carpet': [{'aaxis': {'endlinecolor': '#2a3f5f',
                                                                        'gridcolor': 'white',
                                                                        'linecolor': 'white',
                                                                        'minorgridcolor': 'white',
                                                                        'startlinecolor': '#2a3f5f'},
                                                              'baxis': {'endlinecolor': '#2a3f5f',
                                                                        'gridcolor': 'white',
                                                                        'linecolor': 'white',
                                                                        'minorgridcolor': 'white',
                                                                        'startlinecolor': '#2a3f5f'},
                                                              'type': 'carpet'}], 'table': [
                         {'cells': {'fill': {'color': '#EBF0F8'}, 'line': {'color': 'white'}},
                          'header': {'fill': {'color': '#C8D4E3'}, 'line': {'color': 'white'}},
                          'type': 'table'}], 'barpolar': [{'marker': {
                         'line': {'color': '#E5ECF6', 'width': 0.5},
                         'pattern': {'fillmode': 'overlay', 'size': 10, 'solidity': 0.2}},
                                                           'type': 'barpolar'}],
                                                  'pie': [{'automargin': True, 'type': 'pie'}]},
                                         'layout': {'autotypenumbers': 'strict',
                                                    'colorway': ['#636efa', '#EF553B', '#00cc96',
                                                                 '#ab63fa', '#FFA15A', '#19d3f3',
                                                                 '#FF6692', '#B6E880', '#FF97FF',
                                                                 '#FECB52'],
                                                    'font': {'color': '#2a3f5f'},
                                                    'hovermode': 'closest',
                                                    'hoverlabel': {'align': 'left'},
                                                    'paper_bgcolor': 'white',
                                                    'plot_bgcolor': '#E5ECF6',
                                                    'polar': {'bgcolor': '#E5ECF6',
                                                              'angularaxis': {'gridcolor': 'white',
                                                                              'linecolor': 'white',
                                                                              'ticks': ''},
                                                              'radialaxis': {'gridcolor': 'white',
                                                                             'linecolor': 'white',
                                                                             'ticks': ''}},
                                                    'ternary': {'bgcolor': '#E5ECF6',
                                                                'aaxis': {'gridcolor': 'white',
                                                                          'linecolor': 'white',
                                                                          'ticks': ''},
                                                                'baxis': {'gridcolor': 'white',
                                                                          'linecolor': 'white',
                                                                          'ticks': ''},
                                                                'caxis': {'gridcolor': 'white',
                                                                          'linecolor': 'white',
                                                                          'ticks': ''}},
                                                    'coloraxis': {'colorbar': {'outlinewidth': 0,
                                                                               'ticks': ''}},
                                                    'colorscale': {'sequential': [[0.0, '#0d0887'],
                                                                                  [
                                                                                      0.1111111111111111,
                                                                                      '#46039f'], [
                                                                                      0.2222222222222222,
                                                                                      '#7201a8'], [
                                                                                      0.3333333333333333,
                                                                                      '#9c179e'], [
                                                                                      0.4444444444444444,
                                                                                      '#bd3786'], [
                                                                                      0.5555555555555556,
                                                                                      '#d8576b'], [
                                                                                      0.6666666666666666,
                                                                                      '#ed7953'], [
                                                                                      0.7777777777777778,
                                                                                      '#fb9f3a'], [
                                                                                      0.8888888888888888,
                                                                                      '#fdca26'],
                                                                                  [1.0, '#f0f921']],
                                                                   'sequentialminus': [
                                                                       [0.0, '#0d0887'],
                                                                       [0.1111111111111111,
                                                                        '#46039f'],
                                                                       [0.2222222222222222,
                                                                        '#7201a8'],
                                                                       [0.3333333333333333,
                                                                        '#9c179e'],
                                                                       [0.4444444444444444,
                                                                        '#bd3786'],
                                                                       [0.5555555555555556,
                                                                        '#d8576b'],
                                                                       [0.6666666666666666,
                                                                        '#ed7953'],
                                                                       [0.7777777777777778,
                                                                        '#fb9f3a'],
                                                                       [0.8888888888888888,
                                                                        '#fdca26'],
                                                                       [1.0, '#f0f921']],
                                                                   'diverging': [[0, '#8e0152'],
                                                                                 [0.1, '#c51b7d'],
                                                                                 [0.2, '#de77ae'],
                                                                                 [0.3, '#f1b6da'],
                                                                                 [0.4, '#fde0ef'],
                                                                                 [0.5, '#f7f7f7'],
                                                                                 [0.6, '#e6f5d0'],
                                                                                 [0.7, '#b8e186'],
                                                                                 [0.8, '#7fbc41'],
                                                                                 [0.9, '#4d9221'],
                                                                                 [1, '#276419']]},
                                                    'xaxis': {'gridcolor': 'white',
                                                              'linecolor': 'white', 'ticks': '',
                                                              'title': {'standoff': 15},
                                                              'zerolinecolor': 'white',
                                                              'automargin': True,
                                                              'zerolinewidth': 2},
                                                    'yaxis': {'gridcolor': 'white',
                                                              'linecolor': 'white', 'ticks': '',
                                                              'title': {'standoff': 15},
                                                              'zerolinecolor': 'white',
                                                              'automargin': True,
                                                              'zerolinewidth': 2}, 'scene': {
                                                 'xaxis': {'backgroundcolor': '#E5ECF6',
                                                           'gridcolor': 'white',
                                                           'linecolor': 'white',
                                                           'showbackground': True, 'ticks': '',
                                                           'zerolinecolor': 'white',
                                                           'gridwidth': 2},
                                                 'yaxis': {'backgroundcolor': '#E5ECF6',
                                                           'gridcolor': 'white',
                                                           'linecolor': 'white',
                                                           'showbackground': True, 'ticks': '',
                                                           'zerolinecolor': 'white',
                                                           'gridwidth': 2},
                                                 'zaxis': {'backgroundcolor': '#E5ECF6',
                                                           'gridcolor': 'white',
                                                           'linecolor': 'white',
                                                           'showbackground': True, 'ticks': '',
                                                           'zerolinecolor': 'white',
                                                           'gridwidth': 2}},
                                                    'shapedefaults': {'line': {'color': '#2a3f5f'}},
                                                    'annotationdefaults': {'arrowcolor': '#2a3f5f',
                                                                           'arrowhead': 0,
                                                                           'arrowwidth': 1},
                                                    'geo': {'bgcolor': 'white',
                                                            'landcolor': '#E5ECF6',
                                                            'subunitcolor': 'white',
                                                            'showland': True, 'showlakes': True,
                                                            'lakecolor': 'white'},
                                                    'title': {'x': 0.05},
                                                    'mapbox': {'style': 'light'}}}, 'annotations': [
                     {'font': {'size': 30.0}, 'showarrow': False, 'text': 'Significant p-values',
                      'x': 0.135, 'xanchor': 'center', 'xref': 'paper', 'y': 1.0,
                      'yanchor': 'bottom', 'yref': 'paper'},
                     {'font': {'size': 30.0}, 'showarrow': False,
                      'text': 'out : Classes enrichment representation', 'x': 0.685,
                      'xanchor': 'center', 'xref': 'paper', 'y': 1.0, 'yanchor': 'bottom',
                      'yref': 'paper'}], 'font': {'color': '#111111', 'size': 20},
                            'paper_bgcolor': 'rgba(255, 255, 255, 0)'}}

        self.assertEqual(fig, w_fig)