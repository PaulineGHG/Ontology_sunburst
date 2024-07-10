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

DATA = {'ID': ['FRAMES', 'cdeeg+__FRAMES', 'cdeeg__cdeeg+__FRAMES', 'cdecf__FRAMES', 'gh__FRAMES',
               'eg__cdeeg__cdeeg+__FRAMES', 'eg__FRAMES', 'cde__cdeeg__cdeeg+__FRAMES',
               'cde__cdecf__FRAMES', 'cf__cdecf__FRAMES', 'h__gh__FRAMES',
               'g__eg__cdeeg__cdeeg+__FRAMES', 'g__eg__FRAMES', 'g__gh__FRAMES',
               'f__cf__cdecf__FRAMES', 'e__eg__cdeeg__cdeeg+__FRAMES',
               'e__cde__cdeeg__cdeeg+__FRAMES', 'e__eg__FRAMES', 'e__cde__cdecf__FRAMES',
               'd__cde__cdecf__FRAMES', 'd__cde__cdeeg__cdeeg+__FRAMES', 'c__cde__cdecf__FRAMES',
               'c__cf__cdecf__FRAMES', 'c__cde__cdeeg__cdeeg+__FRAMES', 'ab__FRAMES',
               'b__ab__FRAMES', 'a__ab__FRAMES'],
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
        line = (dico[IDS][i], dico[PARENT][i], dico[LABEL][i], dico[COUNT][i], dico[REF_COUNT][i])
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
                     PARENT: ['salutations'],
                     LABEL: ['bonjour'],
                     COUNT: [2],
                     REF_COUNT: [8]}
        data = add_value_data(data=base_data, m_id='slt', label='salut', value=0.5, base_value=2.3,
                              parent='salutations')
        wanted_data = {IDS: ['bjr', 'slt'],
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
        w_lines = {('cdecf__FRAMES', 'FRAMES', 'cdecf', 3, 18),
                   ('a__ab__FRAMES', 'ab__FRAMES', 'a', 1, 1),
                   ('g__gh__FRAMES', 'gh__FRAMES', 'g', nan, 7),
                   ('g__eg__cdeeg__cdeeg+__FRAMES', 'eg__cdeeg__cdeeg+__FRAMES', 'g', nan, 7),
                   ('d__cde__cdeeg__cdeeg+__FRAMES', 'cde__cdeeg__cdeeg+__FRAMES', 'd', nan, 4),
                   ('b__ab__FRAMES', 'ab__FRAMES', 'b', 2, 2),
                   ('cdeeg+__FRAMES', 'FRAMES', 'cdeeg+', 3, 19),
                   ('e__cde__cdeeg__cdeeg+__FRAMES', 'cde__cdeeg__cdeeg+__FRAMES', 'e', nan, 5),
                   ('e__cde__cdecf__FRAMES', 'cde__cdecf__FRAMES', 'e', nan, 5),
                   ('cde__cdecf__FRAMES', 'cdecf__FRAMES', 'cde', 3, 12),
                   ('c__cde__cdecf__FRAMES', 'cde__cdecf__FRAMES', 'c', 3, 3),
                   ('e__eg__cdeeg__cdeeg+__FRAMES', 'eg__cdeeg__cdeeg+__FRAMES', 'e', nan, 5),
                   ('cf__cdecf__FRAMES', 'cdecf__FRAMES', 'cf', 3, 9),
                   ('ab__FRAMES', 'FRAMES', 'ab', 3, 3),
                   ('g__eg__FRAMES', 'eg__FRAMES', 'g', nan, 7),
                   ('cdeeg__cdeeg+__FRAMES', 'cdeeg+__FRAMES', 'cdeeg', 3, 19),
                   ('h__gh__FRAMES', 'gh__FRAMES', 'h', nan, 8),
                   ('gh__FRAMES', 'FRAMES', 'gh', nan, 15),
                   ('c__cde__cdeeg__cdeeg+__FRAMES', 'cde__cdeeg__cdeeg+__FRAMES', 'c', 3, 3),
                   ('f__cf__cdecf__FRAMES', 'cf__cdecf__FRAMES', 'f', nan, 6),
                   ('e__eg__FRAMES', 'eg__FRAMES', 'e', nan, 5),
                   ('d__cde__cdecf__FRAMES', 'cde__cdecf__FRAMES', 'd', nan, 4),
                   ('eg__FRAMES', 'FRAMES', 'eg', nan, 12),
                   ('c__cf__cdecf__FRAMES', 'cf__cdecf__FRAMES', 'c', 3, 3),
                   ('eg__cdeeg__cdeeg+__FRAMES', 'cdeeg__cdeeg+__FRAMES', 'eg', nan, 12),
                   ('FRAMES', '', 'FRAMES', 6, 36),
                   ('cde__cdeeg__cdeeg+__FRAMES', 'cdeeg__cdeeg+__FRAMES', 'cde', 3, 12)}
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
