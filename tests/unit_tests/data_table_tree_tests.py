import unittest
from unittest.mock import patch
import io
import sys

import numpy as np
import pandas as pd
from functools import wraps
from ontosunburst.data_table_tree import *
from ontosunburst.ontology import *
from ontosunburst.ontosunburst import ontosunburst

"""
Tests manually good file creation.
No automatic tests integrated.
"""

# ==================================================================================================
# GLOBAL
# ==================================================================================================

# GENERAL DICT ONTO (METACYC, KEGG)
# --------------------------------------------------------------------------------------------------

MET_LST = ['a', 'b', 'c']
MET_REF = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
MET_LAB = [1, 2, 3]
MET_RAB = [1, 2, 3, 4, 5, 6, 7, 8]
MET_RAB3 = [8, 7, 6, 5, 4, 3, 2, 1]
MC_ONTO = {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf'], 'd': ['cde'], 'e': ['cde', 'eg'],
           'f': ['cf'], 'g': ['gh', 'eg'], 'h': ['gh'],
           'ab': [ROOTS[METACYC]], 'cde': ['cdecf', 'cdeeg'], 'cf': ['cdecf'],
           'eg': ['cdeeg', ROOTS[METACYC]], 'gh': [ROOTS[METACYC]],
           'cdecf': [ROOTS[METACYC]], 'cdeeg': ['cdeeg+'], 'cdeeg+': [ROOTS[METACYC]]}

MC_ONTO_CH = {'ab': ['a', 'b'], 'cde': ['c', 'd', 'e'], 'cf': ['c', 'f'], 'eg': ['e', 'g'],
              'gh': ['g', 'h'], 'FRAMES': ['ab', 'eg', 'gh', 'cdecf', 'cdeeg+'],
              'cdecf': ['cde', 'cf'], 'cdeeg': ['cde', 'eg'], 'cdeeg+': ['cdeeg']}

MET_LAB_D = {'FRAMES': 6, 'cde': 3, 'cf': 3, 'cdecf': 3, 'cdeeg+': 3, 'cdeeg': 3, 'c': 3, 'ab': 3,
             'b': 2, 'a': 1}
MET_RAB_D = {'FRAMES': 36, 'cdeeg+': 19, 'cdeeg': 19, 'cdecf': 18, 'gh': 15, 'eg': 12, 'cde': 12,
             'cf': 9, 'h': 8, 'g': 7, 'f': 6, 'e': 5, 'd': 4, 'c': 3, 'ab': 3, 'b': 2, 'a': 1}


# ontosunburst(['c', 'd', 'e', 'f', 'cf'], abundances=[3, 4, 5, 2, 2], class_ontology=MC_ONTO,
#              root=ROOTS[METACYC], output='here', show_leaves=True, full=False)
# ontosunburst(metabolic_objects=MET_LST, abundances=MET_LAB, ref_base=True,
#              reference_set=MET_REF, ref_abundances=MET_RAB3, class_ontology=MC_ONTO,
#              root=ROOTS[METACYC], output='here', show_leaves=True, full=False)


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
    df = pd.DataFrame(data=dico.values(), index=dico.keys())
    for l in df.T.iterrows():
        line = (l[1][IDS], l[1][PARENT], l[1][LABEL], l[1][COUNT], l[1][REF_COUNT])
        lines.add(line)
    return lines


# ==================================================================================================
# UNIT TESTS
# ==================================================================================================

# TEST
# --------------------------------------------------------------------------------------------------

class TestDataTable(unittest.TestCase):

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
        np.testing.assert_equal(sub_abu, np.nan)

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
                   ('g__gh__FRAMES', 'gh__FRAMES', 'g', np.nan, 7),
                   ('g__eg__cdeeg__cdeeg+__FRAMES', 'eg__cdeeg__cdeeg+__FRAMES', 'g', np.nan, 7),
                   ('d__cde__cdeeg__cdeeg+__FRAMES', 'cde__cdeeg__cdeeg+__FRAMES', 'd', np.nan, 4),
                   ('b__ab__FRAMES', 'ab__FRAMES', 'b', 2, 2),
                   ('cdeeg+__FRAMES', 'FRAMES', 'cdeeg+', 3, 19),
                   ('e__cde__cdeeg__cdeeg+__FRAMES', 'cde__cdeeg__cdeeg+__FRAMES', 'e', np.nan, 5),
                   ('e__cde__cdecf__FRAMES', 'cde__cdecf__FRAMES', 'e', np.nan, 5),
                   ('cde__cdecf__FRAMES', 'cdecf__FRAMES', 'cde', 3, 12),
                   ('c__cde__cdecf__FRAMES', 'cde__cdecf__FRAMES', 'c', 3, 3),
                   ('e__eg__cdeeg__cdeeg+__FRAMES', 'eg__cdeeg__cdeeg+__FRAMES', 'e', np.nan, 5),
                   ('cf__cdecf__FRAMES', 'cdecf__FRAMES', 'cf', 3, 9),
                   ('ab__FRAMES', 'FRAMES', 'ab', 3, 3),
                   ('g__eg__FRAMES', 'eg__FRAMES', 'g', np.nan, 7),
                   ('cdeeg__cdeeg+__FRAMES', 'cdeeg+__FRAMES', 'cdeeg', 3, 19),
                   ('h__gh__FRAMES', 'gh__FRAMES', 'h', np.nan, 8),
                   ('gh__FRAMES', 'FRAMES', 'gh', np.nan, 15),
                   ('c__cde__cdeeg__cdeeg+__FRAMES', 'cde__cdeeg__cdeeg+__FRAMES', 'c', 3, 3),
                   ('f__cf__cdecf__FRAMES', 'cf__cdecf__FRAMES', 'f', np.nan, 6),
                   ('e__eg__FRAMES', 'eg__FRAMES', 'e', np.nan, 5),
                   ('d__cde__cdecf__FRAMES', 'cde__cdecf__FRAMES', 'd', np.nan, 4),
                   ('eg__FRAMES', 'FRAMES', 'eg', np.nan, 12),
                   ('c__cf__cdecf__FRAMES', 'cf__cdecf__FRAMES', 'c', 3, 3),
                   ('eg__cdeeg__cdeeg+__FRAMES', 'cdeeg__cdeeg+__FRAMES', 'eg', np.nan, 12),
                   ('FRAMES', '', 'FRAMES', 6, 36),
                   ('cde__cdeeg__cdeeg+__FRAMES', 'cdeeg__cdeeg+__FRAMES', 'cde', 3, 12)}
        self.assertEqual(lines, w_lines)
