import copy
import json
import os.path
import unittest
import io
from functools import wraps

from numpy import nan
from ontosunburst.sunburst_fig import *

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
ENRICH_P_VAL = [0.0, -5.420266413988895, 2.509702991379166, 2.948803113091024, -1.2341542222355069,
                -1.103304935668835, nan, nan, 0.4034095751193356, 0.0]

ENRICH_REF_AB = {'00': 100, '01': 40, '02': 30, '03': 20, '04': 10, '05': 20, '06': 5, '07': 1,
                 '08': 1, '09': 3}

DATA = {'ID': ['FRAMES', 'cdeeg+__FRAMES', 'cdeeg__cdeeg+__FRAMES', 'cdecf__FRAMES', 'gh__FRAMES',
               'eg__cdeeg__cdeeg+__FRAMES', 'eg__FRAMES', 'cde__cdeeg__cdeeg+__FRAMES',
               'cde__cdecf__FRAMES', 'cf__cdecf__FRAMES', 'h__gh__FRAMES',
               'g__eg__cdeeg__cdeeg+__FRAMES', 'g__eg__FRAMES', 'g__gh__FRAMES',
               'f__cf__cdecf__FRAMES', 'e__eg__cdeeg__cdeeg+__FRAMES',
               'e__cde__cdeeg__cdeeg+__FRAMES', 'e__eg__FRAMES', 'e__cde__cdecf__FRAMES',
               'd__cde__cdecf__FRAMES', 'd__cde__cdeeg__cdeeg+__FRAMES', 'c__cde__cdecf__FRAMES',
               'c__cf__cdecf__FRAMES', 'c__cde__cdeeg__cdeeg+__FRAMES', 'ab__FRAMES',
               'b__ab__FRAMES', 'a__ab__FRAMES'],
        'Onto ID': ['FRAMES', 'cdeeg+', 'cdeeg', 'cdecf', 'gh', 'eg', 'eg', 'cde', 'cde', 'cf',
                    'h', 'g', 'g', 'g', 'f', 'e', 'e', 'e', 'e', 'd', 'd', 'c', 'c', 'c', 'ab',
                    'b', 'a'],
        'Parent': ['', 'FRAMES', 'cdeeg+__FRAMES', 'FRAMES', 'FRAMES', 'cdeeg__cdeeg+__FRAMES',
                   'FRAMES', 'cdeeg__cdeeg+__FRAMES', 'cdecf__FRAMES', 'cdecf__FRAMES',
                   'gh__FRAMES', 'eg__cdeeg__cdeeg+__FRAMES', 'eg__FRAMES', 'gh__FRAMES',
                   'cf__cdecf__FRAMES', 'eg__cdeeg__cdeeg+__FRAMES', 'cde__cdeeg__cdeeg+__FRAMES',
                   'eg__FRAMES', 'cde__cdecf__FRAMES', 'cde__cdecf__FRAMES',
                   'cde__cdeeg__cdeeg+__FRAMES', 'cde__cdecf__FRAMES', 'cf__cdecf__FRAMES',
                   'cde__cdeeg__cdeeg+__FRAMES', 'FRAMES', 'ab__FRAMES', 'ab__FRAMES'],
        'Label': ['FRAMES', 'cdeeg+', 'cdeeg', 'cdecf', 'gh', 'eg', 'eg', 'cde', 'cde', 'cf', 'h',
                  'g', 'g', 'g', 'f', 'e', 'e', 'e', 'e', 'd', 'd', 'c', 'c', 'c', 'ab', 'b', 'a'],
        'Count': [6, 3, 3, 3, nan, nan, nan, 3, 3, 3, nan, nan, nan, nan, nan, nan, nan, nan, nan,
                  nan, nan, 3, 3, 3, 3, 2, 1],
        'Reference count': [36, 19, 19, 18, 15, 12, 12, 12, 12, 9, 8, 7, 7, 7, 6, 5, 5, 5, 5, 4, 4,
                            3, 3, 3, 3, 2, 1],
        'Proportion': [1.0, 0.5, 0.5, 0.5, nan, nan, nan, 0.5, 0.5, 0.5, nan, nan, nan, nan, nan,
                       nan, nan, nan, nan, nan, nan, 0.5, 0.5, 0.5, 0.5, 0.3333333333333333,
                       0.16666666666666666],
        'Reference proportion': [1.0, 0.5277777777777778, 0.5277777777777778, 0.5,
                                 0.4166666666666667, 0.3333333333333333, 0.3333333333333333,
                                 0.3333333333333333, 0.3333333333333333, 0.25, 0.2222222222222222,
                                 0.19444444444444445, 0.19444444444444445, 0.19444444444444445,
                                 0.16666666666666666, 0.1388888888888889, 0.1388888888888889,
                                 0.1388888888888889, 0.1388888888888889, 0.1111111111111111,
                                 0.1111111111111111, 0.08333333333333333, 0.08333333333333333,
                                 0.08333333333333333, 0.08333333333333333, 0.05555555555555555,
                                 0.027777777777777776],
        'Relative proportion': [1000000, 283582, 283582, 268656, 223880, 141791, 179104,
                                141791, 153517, 115138, 119402, 82711, 104477, 104477, 76758,
                                59079, 59079, 74626, 63965, 51172, 47263, 38379, 38379, 35447,
                                44776, 29850, 14925]}


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


def save_fig_json(fig, file):
    fig = fig.to_dict()
    with open(file, 'w') as f:
        json.dump(fig, f)


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
        data = copy.deepcopy(ENRICH_DATA)
        data[PVAL] = ENRICH_P_VAL
        text_list = get_hover_fig_text(data, ENRICHMENT_A, True)
        self.assertEqual(len(text_list), 10)
        self.assertEqual(text_list[0], 'P value: 1.0<br>Count: <b>50</b>'
                                       '<br>Reference count: 100<br>Proportion: <b>100%</b>'
                                       '<br>Reference proportion: 100%<br>ID: 00')
        self.assertEqual(text_list[5], 'P value: 0.07883064215278136<br>Count: <b>5</b>'
                                       '<br>Reference count: 20<br>Proportion: <b>10.0%</b>'
                                       '<br>Reference proportion: 20.0%<br>ID: 05')

    @test_for(get_hover_fig_text)
    def test_get_hover_fig_text_enrich_no_ref(self):
        data = copy.deepcopy(ENRICH_DATA)
        data[PVAL] = ENRICH_P_VAL
        text_list = get_hover_fig_text(data, ENRICHMENT_A, False)
        self.assertEqual(len(text_list), 10)
        self.assertEqual(text_list[0], 'P value: 1.0<br>Count: <b>50</b>'
                                       '<br>Proportion: <b>100%</b>'
                                       '<br>ID: 00')
        self.assertEqual(text_list[5], 'P value: 0.07883064215278136<br>Count: <b>5</b>'
                                       '<br>Proportion: <b>10.0%</b>'
                                       '<br>ID: 05')

    @test_for(get_hover_fig_text)
    def test_get_hover_fig_text_topology_ref(self):
        data = copy.deepcopy(ENRICH_DATA)
        text_list = get_hover_fig_text(data, TOPOLOGY_A, True)
        self.assertEqual(len(text_list), 10)
        self.assertEqual(text_list[0], 'Count: <b>50</b><br>Reference count: 100'
                                       '<br>Proportion: <b>100%</b>'
                                       '<br>Reference proportion: 100%<br>ID: 00')
        self.assertEqual(text_list[5], 'Count: <b>5</b><br>Reference count: 20'
                                       '<br>Proportion: <b>10.0%</b>'
                                       '<br>Reference proportion: 20.0%<br>ID: 05')

    @test_for(get_hover_fig_text)
    def test_get_hover_fig_text_topology_no_ref(self):
        data = copy.deepcopy(ENRICH_DATA)
        text_list = get_hover_fig_text(data, TOPOLOGY_A, False)
        self.assertEqual(len(text_list), 10)
        self.assertEqual(text_list[0], 'Count: <b>50</b><br>Proportion: <b>100%</b><br>ID: 00')
        self.assertEqual(text_list[5], 'Count: <b>5</b><br>Proportion: <b>10.0%</b><br>ID: 05')

    @test_for(generate_sunburst_fig)
    def test_generate_sunburst_fig_case1(self):
        data = copy.deepcopy(ENRICH_DATA)
        fig = generate_sunburst_fig(data, 'case1', analysis=ENRICHMENT_A, write_fig=False,
                                    ref_classes_abundance=ENRICH_REF_AB, test=HYPERGEO_TEST)
        w_fig_file = os.path.join('test_files', 'fig_case1.json')
        fig = json.dumps(fig.to_dict(), sort_keys=True)
        with open(w_fig_file, 'r') as f:
            w_fig = json.dumps(json.load(f), sort_keys=True)
        self.assertEqual(fig, w_fig)

    @test_for(generate_sunburst_fig)
    def test_generate_sunburst_fig_case2(self):
        data = copy.deepcopy(ENRICH_DATA)
        fig = generate_sunburst_fig(data, 'case2', analysis=ENRICHMENT_A, write_fig=False,
                                    ref_classes_abundance=ENRICH_REF_AB, test=BINOMIAL_TEST)
        w_fig_file = os.path.join('test_files', 'fig_case2.json')
        fig = json.dumps(fig.to_dict(), sort_keys=True)
        with open(w_fig_file, 'r') as f:
            w_fig = json.dumps(json.load(f), sort_keys=True)
        self.assertEqual(fig, w_fig)

    @test_for(generate_sunburst_fig)
    def test_generate_sunburst_fig_case3(self):
        data = copy.deepcopy(ENRICH_DATA)
        fig = generate_sunburst_fig(data, 'case3', analysis=ENRICHMENT_A,
                                    ref_classes_abundance=ENRICH_REF_AB, test=HYPERGEO_TEST,
                                    root_cut=ROOT_UNCUT, write_fig=False,
                                    title='Another title', colorscale='PuOr_r',
                                    bg_color='#222222', font_color='#eeeeee', font_size=25,
                                    table_legend='Number')
        w_fig_file = os.path.join('test_files', 'fig_case3.json')
        fig = json.dumps(fig.to_dict(), sort_keys=True)
        with open(w_fig_file, 'r') as f:
            w_fig = json.dumps(json.load(f), sort_keys=True)
        self.assertEqual(fig, w_fig)

    @test_for(generate_sunburst_fig)
    def test_generate_sunburst_fig_case4(self):
        data = copy.deepcopy(DATA)
        fig = generate_sunburst_fig(data, 'case4', analysis=TOPOLOGY_A, bg_color='black',
                                    font_color='white', write_fig=False)
        w_fig_file = os.path.join('test_files', 'fig_case4.json')
        fig = json.dumps(fig.to_dict(), sort_keys=True)
        with open(w_fig_file, 'r') as f:
            w_fig = json.dumps(json.load(f), sort_keys=True)
        self.assertEqual(fig, w_fig)



