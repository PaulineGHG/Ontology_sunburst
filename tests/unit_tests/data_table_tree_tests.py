import unittest
import io

from numpy import nan
from functools import wraps
from ontosunburst.data_table_tree import *
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

MC_ONTO = {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf'], 'd': ['cde'], 'e': ['cde', 'eg'],
           'f': ['cf'], 'g': ['gh', 'eg'], 'h': ['gh'],
           'ab': [ROOTS[METACYC]], 'cde': ['cdecf', 'cdeeg'], 'cf': ['cdecf'],
           'eg': ['cdeeg', ROOTS[METACYC]], 'gh': [ROOTS[METACYC]],
           'cdecf': [ROOTS[METACYC]], 'cdeeg': ['cdeeg+'], 'cdeeg+': [ROOTS[METACYC]]}

MET_LAB_D = {'FRAMES': 6, 'cde': 3, 'cf': 3, 'cdecf': 3, 'cdeeg+': 3, 'cdeeg': 3, 'c': 3, 'ab': 3,
             'b': 2, 'a': 1}
MET_RAB_D = {'FRAMES': 36, 'cdeeg+': 19, 'cdeeg': 19, 'cdecf': 18, 'gh': 15, 'eg': 12, 'cde': 12,
             'cf': 9, 'h': 8, 'g': 7, 'f': 6, 'e': 5, 'd': 4, 'c': 3, 'ab': 3, 'b': 2, 'a': 1}

NAMES = {'FRAMES': 'Root', 'cdeeg+': 'CDEEG+', 'cdeeg': 'CDEEG', 'cdecf': 'CDECF', 'gh': 'GH',
         'eg': 'EG', 'cde': 'CDE', 'cf': 'CF', 'h': 'H', 'g': 'G', 'f': 'F', 'e': 'E', 'd': 'D',
         'c': 'C', 'ab': 'AB', 'b': 'B'}

DATA = {'ID': ['FRAMES', 'cdeeg+__FRAMES', 'cdeeg__cdeeg+__FRAMES', 'cdecf__FRAMES', 'gh__FRAMES',
               'eg__cdeeg__cdeeg+__FRAMES', 'eg__FRAMES', 'cde__cdeeg__cdeeg+__FRAMES',
               'cde__cdecf__FRAMES', 'cf__cdecf__FRAMES', 'h__gh__FRAMES',
               'g__eg__cdeeg__cdeeg+__FRAMES', 'g__eg__FRAMES', 'g__gh__FRAMES',
               'f__cf__cdecf__FRAMES', 'e__eg__cdeeg__cdeeg+__FRAMES',
               'e__cde__cdeeg__cdeeg+__FRAMES', 'e__eg__FRAMES', 'e__cde__cdecf__FRAMES',
               'd__cde__cdecf__FRAMES', 'd__cde__cdeeg__cdeeg+__FRAMES', 'c__cde__cdecf__FRAMES',
               'c__cf__cdecf__FRAMES', 'c__cde__cdeeg__cdeeg+__FRAMES', 'ab__FRAMES',
               'b__ab__FRAMES', 'a__ab__FRAMES'],
        'Onto ID': ['FRAMES', 'cdeeg+', 'cdeeg', 'cdecf', 'gh', 'eg', 'eg', 'cde', 'cde', 'cf', 'h',
                    'g', 'g', 'g', 'f', 'e', 'e', 'e', 'e', 'd', 'd', 'c', 'c', 'c', 'ab', 'b', 'a'],
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
                            3, 3, 3, 3, 2, 1]}

W_PROP = [1.0, 0.5, 0.5, 0.5, nan, nan, nan, 0.5, 0.5, 0.5, nan, nan, nan, nan, nan, nan, nan, nan,
          nan, nan, nan, 0.5, 0.5, 0.5, 0.5, 0.3333333333333333, 0.16666666666666666]

W_REF_PROP = [1.0, 0.5277777777777778, 0.5277777777777778, 0.5, 0.4166666666666667,
              0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.25,
              0.2222222222222222, 0.19444444444444445, 0.19444444444444445, 0.19444444444444445,
              0.16666666666666666, 0.1388888888888889, 0.1388888888888889, 0.1388888888888889,
              0.1388888888888889, 0.1111111111111111, 0.1111111111111111, 0.08333333333333333,
              0.08333333333333333, 0.08333333333333333, 0.08333333333333333, 0.05555555555555555,
              0.027777777777777776]

W_RELAT_PROP = [1000000, 283582, 283582, 268656, 223880, 141791, 179104, 141791, 153517, 115138,
                119402, 82711, 104477, 104477, 76758, 59079, 59079, 74626, 63965, 51172, 47263,
                38379, 38379, 35447, 44776, 29850, 14925]

DATA_PROP = DATA
DATA_PROP[PROP] = W_PROP
DATA_PROP[REF_PROP] = W_REF_PROP

DATA_R_PROP = DATA_PROP
DATA_R_PROP[RELAT_PROP] = W_RELAT_PROP


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


def data_to_lines(dico):
    lines = set()
    for i in range(len(dico[IDS])):
        line = (dico[IDS][i], dico[ONTO_ID][i], dico[PARENT][i], dico[LABEL][i], dico[COUNT][i],
                dico[REF_COUNT][i])
        if PROP in dico:
            line = line + (dico[PROP][i],)
        if REF_PROP in dico:
            line = line + (dico[REF_PROP][i],)
        if RELAT_PROP in dico:
            line = line + (dico[RELAT_PROP][i],)
        if PVAL in dico:
            line = line + (dico[PVAL][i],)
        lines.add(line)
    return lines


# ==================================================================================================
# UNIT TESTS
# ==================================================================================================

# TEST
# --------------------------------------------------------------------------------------------------

class TestGenerateDataTable(unittest.TestCase):

    @test_for(get_sub_abundance)
    def test_get_sub_abundances_exists_diff(self):
        sub_abu = get_sub_abundance(MET_LAB_D, 'cf', 9)
        self.assertEqual(sub_abu, 3)

    @test_for(get_sub_abundance)
    def test_get_sub_abundances_exists_equ(self):
        sub_abu = get_sub_abundance(MET_LAB_D, 'a', 1)
        self.assertEqual(sub_abu, 1)

    @test_for(get_sub_abundance)
    def test_get_sub_abundances_not_exists(self):
        sub_abu = get_sub_abundance(MET_LAB_D, 'eg', 12)
        self.assertTrue(np.isnan(sub_abu))

    @test_for(get_sub_abundance)
    def test_get_sub_abundances_no_sub(self):
        sub_abu = get_sub_abundance(None, 'eg', 12)
        self.assertEqual(sub_abu, 12)

    @test_for(add_value_data)
    def test_add_value_data(self):
        base_data = {IDS: ['bjr'],
                     ONTO_ID: ['Bjr_0'],
                     PARENT: ['salutations'],
                     LABEL: ['bonjour'],
                     COUNT: [2],
                     REF_COUNT: [8]}
        data = add_value_data(data=base_data, m_id='slt', onto_id='sl_1', label='salut', value=0.5,
                              base_value=2.3, parent='salutations')
        wanted_data = {IDS: ['bjr', 'slt'],
                       ONTO_ID: ['Bjr_0', 'sl_1'],
                       PARENT: ['salutations', 'salutations'],
                       LABEL: ['bonjour', 'salut'],
                       COUNT: [2, 0.5],
                       REF_COUNT: [8, 2.3]}
        self.assertEqual(data, wanted_data)

    @test_for(get_all_ids)
    def test_get_all_c_ids(self):
        all_ids = get_all_ids('c', 'c', MC_ONTO, ROOTS[METACYC], set())
        wanted_ids = {'c__cf__cdecf__FRAMES', 'c__cde__cdecf__FRAMES',
                      'c__cde__cdeeg__cdeeg+__FRAMES'}
        self.assertEqual(all_ids, wanted_ids)

    @test_for(get_all_ids)
    def test_get_all_e_ids(self):
        all_ids = get_all_ids('e', 'e', MC_ONTO, ROOTS[METACYC], set())
        wanted_ids = {'e__cde__cdeeg__cdeeg+__FRAMES', 'e__eg__FRAMES',
                      'e__eg__cdeeg__cdeeg+__FRAMES', 'e__cde__cdecf__FRAMES'}
        self.assertEqual(all_ids, wanted_ids)

    @test_for(get_all_ids)
    def test_get_all_eg_ids(self):
        all_ids = get_all_ids('eg', 'eg', MC_ONTO, ROOTS[METACYC], set())
        wanted_ids = {'eg__FRAMES', 'eg__cdeeg__cdeeg+__FRAMES'}
        self.assertEqual(all_ids, wanted_ids)

    @test_for(get_fig_parameters)
    def test_get_fig_parameters(self):
        data = get_fig_parameters(classes_abondance=MET_RAB_D, parent_dict=MC_ONTO,
                                  root_item=ROOTS[METACYC], subset_abundance=MET_LAB_D, names=None)
        lines = data_to_lines(data)
        w_lines = {('c__cf__cdecf__FRAMES', 'c', 'cf__cdecf__FRAMES', 'c', 3, 3),
                   ('b__ab__FRAMES', 'b', 'ab__FRAMES', 'b', 2, 2),
                   ('eg__FRAMES', 'eg', 'FRAMES', 'eg', nan, 12),
                   ('g__eg__cdeeg__cdeeg+__FRAMES', 'g', 'eg__cdeeg__cdeeg+__FRAMES', 'g', nan, 7),
                   ('h__gh__FRAMES', 'h', 'gh__FRAMES', 'h', nan, 8),
                   ('f__cf__cdecf__FRAMES', 'f', 'cf__cdecf__FRAMES', 'f', nan, 6),
                   ('c__cde__cdecf__FRAMES', 'c', 'cde__cdecf__FRAMES', 'c', 3, 3),
                   ('e__cde__cdeeg__cdeeg+__FRAMES', 'e', 'cde__cdeeg__cdeeg+__FRAMES', 'e', nan, 5),
                   ('gh__FRAMES', 'gh', 'FRAMES', 'gh', nan, 15),
                   ('e__eg__cdeeg__cdeeg+__FRAMES', 'e', 'eg__cdeeg__cdeeg+__FRAMES', 'e', nan, 5),
                   ('g__gh__FRAMES', 'g', 'gh__FRAMES', 'g', nan, 7),
                   ('d__cde__cdecf__FRAMES', 'd', 'cde__cdecf__FRAMES', 'd', nan, 4),
                   ('cdecf__FRAMES', 'cdecf', 'FRAMES', 'cdecf', 3, 18),
                   ('cf__cdecf__FRAMES', 'cf', 'cdecf__FRAMES', 'cf', 3, 9),
                   ('g__eg__FRAMES', 'g', 'eg__FRAMES', 'g', nan, 7),
                   ('ab__FRAMES', 'ab', 'FRAMES', 'ab', 3, 3),
                   ('e__cde__cdecf__FRAMES', 'e', 'cde__cdecf__FRAMES', 'e', nan, 5),
                   ('d__cde__cdeeg__cdeeg+__FRAMES', 'd', 'cde__cdeeg__cdeeg+__FRAMES', 'd', nan, 4),
                   ('eg__cdeeg__cdeeg+__FRAMES', 'eg', 'cdeeg__cdeeg+__FRAMES', 'eg', nan, 12),
                   ('a__ab__FRAMES', 'a', 'ab__FRAMES', 'a', 1, 1),
                   ('cde__cdeeg__cdeeg+__FRAMES', 'cde', 'cdeeg__cdeeg+__FRAMES', 'cde', 3, 12),
                   ('cdeeg+__FRAMES', 'cdeeg+', 'FRAMES', 'cdeeg+', 3, 19),
                   ('cde__cdecf__FRAMES', 'cde', 'cdecf__FRAMES', 'cde', 3, 12),
                   ('cdeeg__cdeeg+__FRAMES', 'cdeeg', 'cdeeg+__FRAMES', 'cdeeg', 3, 19),
                   ('c__cde__cdeeg__cdeeg+__FRAMES', 'c', 'cde__cdeeg__cdeeg+__FRAMES', 'c', 3, 3),
                   ('e__eg__FRAMES', 'e', 'eg__FRAMES', 'e', nan, 5),
                   ('FRAMES', 'FRAMES', '', 'FRAMES', 6, 36)}
        self.assertEqual(lines, w_lines)

    @test_for(get_fig_parameters)
    def test_get_fig_parameters_names(self):
        data = get_fig_parameters(classes_abondance=MET_RAB_D, parent_dict=MC_ONTO,
                                  root_item=ROOTS[METACYC], subset_abundance=MET_LAB_D, names=NAMES)
        lines = data_to_lines(data)
        w_lines = {('cf__cdecf__FRAMES', 'cf', 'cdecf__FRAMES', 'CF', 3, 9),
                   ('e__eg__FRAMES', 'e', 'eg__FRAMES', 'E', nan, 5),
                   ('cde__cdeeg__cdeeg+__FRAMES', 'cde', 'cdeeg__cdeeg+__FRAMES', 'CDE', 3, 12),
                   ('gh__FRAMES', 'gh', 'FRAMES', 'GH', nan, 15),
                   ('a__ab__FRAMES', 'a', 'ab__FRAMES', 'a', 1, 1),
                   ('g__gh__FRAMES', 'g', 'gh__FRAMES', 'G', nan, 7),
                   ('eg__cdeeg__cdeeg+__FRAMES', 'eg', 'cdeeg__cdeeg+__FRAMES', 'EG', nan, 12),
                   ('e__eg__cdeeg__cdeeg+__FRAMES', 'e', 'eg__cdeeg__cdeeg+__FRAMES', 'E', nan, 5),
                   ('g__eg__FRAMES', 'g', 'eg__FRAMES', 'G', nan, 7),
                   ('FRAMES', 'FRAMES', '', 'FRAMES', 6, 36),
                   ('cdeeg+__FRAMES', 'cdeeg+', 'FRAMES', 'CDEEG+', 3, 19),
                   ('cdeeg__cdeeg+__FRAMES', 'cdeeg', 'cdeeg+__FRAMES', 'CDEEG', 3, 19),
                   ('d__cde__cdecf__FRAMES', 'd', 'cde__cdecf__FRAMES', 'D', nan, 4),
                   ('eg__FRAMES', 'eg', 'FRAMES', 'EG', nan, 12),
                   ('cdecf__FRAMES', 'cdecf', 'FRAMES', 'CDECF', 3, 18),
                   ('e__cde__cdeeg__cdeeg+__FRAMES', 'e', 'cde__cdeeg__cdeeg+__FRAMES', 'E', nan, 5),
                   ('c__cf__cdecf__FRAMES', 'c', 'cf__cdecf__FRAMES', 'C', 3, 3),
                   ('cde__cdecf__FRAMES', 'cde', 'cdecf__FRAMES', 'CDE', 3, 12),
                   ('d__cde__cdeeg__cdeeg+__FRAMES', 'd', 'cde__cdeeg__cdeeg+__FRAMES', 'D', nan, 4),
                   ('c__cde__cdecf__FRAMES', 'c', 'cde__cdecf__FRAMES', 'C', 3, 3),
                   ('f__cf__cdecf__FRAMES', 'f', 'cf__cdecf__FRAMES', 'F', nan, 6),
                   ('g__eg__cdeeg__cdeeg+__FRAMES', 'g', 'eg__cdeeg__cdeeg+__FRAMES', 'G', nan, 7),
                   ('ab__FRAMES', 'ab', 'FRAMES', 'AB', 3, 3),
                   ('c__cde__cdeeg__cdeeg+__FRAMES', 'c', 'cde__cdeeg__cdeeg+__FRAMES', 'C', 3, 3),
                   ('b__ab__FRAMES', 'b', 'ab__FRAMES', 'B', 2, 2),
                   ('h__gh__FRAMES', 'h', 'gh__FRAMES', 'H', nan, 8),
                   ('e__cde__cdecf__FRAMES', 'e', 'cde__cdecf__FRAMES', 'E', nan, 5)}
        self.assertEqual(lines, w_lines)


class TestAddProportionDataTable(unittest.TestCase):

    @test_for(get_data_proportion)
    def test_get_data_proportion_no_relative(self):
        data = get_data_proportion(DATA, False)
        for i in range(len(data[PROP])):
            if np.isnan(data[PROP][i]):
                self.assertTrue(np.isnan(W_PROP[i]))
            else:
                self.assertEqual(data[PROP][i], W_PROP[i])

    @test_for(get_data_proportion)
    def test_get_data_proportion_no_relative_ref(self):
        data = get_data_proportion(DATA, False)
        for i in range(len(data[REF_PROP])):
            self.assertEqual(data[REF_PROP][i], W_REF_PROP[i])

    @test_for(get_data_proportion)
    def test_get_data_proportion_relative(self):
        data = get_data_proportion(DATA, True)
        for i in range(len(data[RELAT_PROP])):
            self.assertEqual(data[RELAT_PROP][i], W_RELAT_PROP[i])

    @test_for(get_relative_prop)
    def test_get_relative_prop(self):
        data = DATA_PROP
        data[RELAT_PROP] = [x for x in data[PROP]]
        data = get_relative_prop(data, '')
        for i in range(len(data[RELAT_PROP])):
            self.assertEqual(data[RELAT_PROP][i], W_RELAT_PROP[i])


ENRICH_DATA = {IDS: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
               ONTO_ID: ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09'],
               PARENT: ['', 0, 0, 0, 0, 1, 1, 1, 2, 2],
               LABEL: ['00', '01', '02', '03', '04', '05', '06', '07', '08', '09'],
               COUNT: [50, 5, 25, 20, 1, 5, nan, nan, 1, 1],
               REF_COUNT: [100, 40, 30, 20, 10, 20, 5, 1, 1, 3],
               PROP: [1, 0.1, 0.5, 0.4, 0.02, 0.1, nan, nan, 0.02, 0.02],
               REF_PROP: [1, 0.4, 0.3, 0.2, 0.1, 0.2, 0.05, 0.01, 0.01, 0.03],
               RELAT_PROP: [1, 0.4, 0.3, 0.2, 0.1, 0.2, 0.05, 0.01, 0.01, 0.03]}
ENRICH_REF_AB = {'00': 100, '01': 40, '02': 30, '03': 20, '04': 10, '05': 20, '06': 5, '07': 1,
                 '08': 1, '09': 3}
LABEL_NAMES = ['r', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine']

# Expected :
# Over : 2, 3 | Under : 1, 4, 5 | No diff : 0, 8, 9 | Nan : 6, 7


class TestEnrichmentAnalysis(unittest.TestCase):

    @test_for(get_data_enrichment_analysis)
    def test_get_data_enrichment_analysis_single_value(self):
        data, significant = get_data_enrichment_analysis(ENRICH_DATA, ENRICH_REF_AB, BINOMIAL_TEST)
        lines = data_to_lines(data)
        p_value_1 = [l[-1] for l in lines if l[0] == 1][0]
        M = 100
        N = 50
        m = 40
        n = 5
        exp_p_value_1 = stats.binomtest(n, N, m / M, alternative='two-sided').pvalue
        exp_p_value_1 = np.log10(exp_p_value_1)
        self.assertEqual(p_value_1, exp_p_value_1)

    @test_for(get_data_enrichment_analysis)
    def test_get_data_enrichment_analysis_binomial(self):
        data, significant = get_data_enrichment_analysis(ENRICH_DATA, ENRICH_REF_AB, BINOMIAL_TEST)
        lines = data_to_lines(data)
        exp_significant = {'01': 3.7996e-06, '03': 0.0011251149, '02': 0.0030924096}
        exp_lines = {(7, '07', 1, '07', nan, 1, nan, 0.01, 0.01, nan),
                     (3, '03', 0, '03', 20, 20, 0.4, 0.2, 0.2, 2.948803113091024),
                     (5, '05', 1, '05', 5, 20, 0.1, 0.2, 0.2, -1.103304935668835),
                     (4, '04', 0, '04', 1, 10, 0.02, 0.1, 0.1, -1.2341542222355069),
                     (2, '02', 0, '02', 25, 30, 0.5, 0.3, 0.3, 2.509702991379166),
                     (1, '01', 0, '01', 5, 40, 0.1, 0.4, 0.4, -5.420266413988895),
                     (9, '09', 2, '09', 1, 3, 0.02, 0.03, 0.03, 0.0),
                     (8, '08', 2, '08', 1, 1, 0.02, 0.01, 0.01, 0.4034095751193356),
                     (0, '00', '', '00', 50, 100, 1, 1, 1, 0.0),
                     (6, '06', 1, '06', nan, 5, nan, 0.05, 0.05, nan)}
        self.assertEqual(lines, exp_lines)
        self.assertEqual(significant, exp_significant)

    @test_for(get_data_enrichment_analysis)
    def test_get_data_enrichment_analysis_hypergeometric(self):
        data, significant = get_data_enrichment_analysis(ENRICH_DATA, ENRICH_REF_AB, HYPERGEO_TEST)
        lines = data_to_lines(data)
        exp_lines = {(8, '08', 2, '08', 1, 1, 0.02, 0.01, 0.01, -0.0),
                     (3, '03', 0, '03', 20, 20, 0.4, 0.2, 0.2, 6.754831139005899),
                     (6, '06', 1, '06', nan, 5, nan, 0.05, 0.05, nan),
                     (5, '05', 1, '05', 5, 20, 0.1, 0.2, 0.2, -1.6413993451973743),
                     (9, '09', 2, '09', 1, 3, 0.02, 0.03, 0.03, -1.4464911998299308e-16),
                     (4, '04', 0, '04', 1, 10, 0.02, 0.1, 0.1, -1.8051946563380086),
                     (2, '02', 0, '02', 25, 30, 0.5, 0.3, 0.3, 4.692610428021241),
                     (0, '00', '', '00', 50, 100, 1, 1, 1, 0.3010299956639812),
                     (1, '01', 0, '01', 5, 40, 0.1, 0.4, 0.4, -9.138873998573988),
                     (7, '07', 1, '07', nan, 1, nan, 0.01, 0.01, nan)}
        exp_significant = {'01': 7e-10, '03': 1.759e-07, '02': 2.0295e-05}
        self.assertEqual(lines, exp_lines)
        self.assertEqual(significant, exp_significant)

    def test_get_data_enrichment_analysis_names(self):
        data = ENRICH_DATA
        data[LABEL] = LABEL_NAMES
        data, significant = get_data_enrichment_analysis(data, ENRICH_REF_AB, BINOMIAL_TEST)
        lines = data_to_lines(data)
        exp_lines = {(2, '02', 0, 'two', 25, 30, 0.5, 0.3, 0.3, 2.509702991379166),
                     (4, '04', 0, 'four', 1, 10, 0.02, 0.1, 0.1, -1.2341542222355069),
                     (9, '09', 2, 'nine', 1, 3, 0.02, 0.03, 0.03, 0.0),
                     (0, '00', '', 'r', 50, 100, 1, 1, 1, 0.0),
                     (6, '06', 1, 'six', nan, 5, nan, 0.05, 0.05, nan),
                     (7, '07', 1, 'seven', nan, 1, nan, 0.01, 0.01, nan),
                     (1, '01', 0, 'one', 5, 40, 0.1, 0.4, 0.4, -5.420266413988895),
                     (5, '05', 1, 'five', 5, 20, 0.1, 0.2, 0.2, -1.103304935668835),
                     (3, '03', 0, 'three', 20, 20, 0.4, 0.2, 0.2, 2.948803113091024),
                     (8, '08', 2, 'eight', 1, 1, 0.02, 0.01, 0.01, 0.4034095751193356)}
        exp_significant = {'01': 3.7996e-06, '03': 0.0011251149, '02': 0.0030924096}
        self.assertEqual(lines, exp_lines)
        self.assertEqual(significant, exp_significant)
