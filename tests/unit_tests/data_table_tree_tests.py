import unittest
from unittest.mock import patch
import io
import sys

import numpy as np
import pandas as pd
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

# MET_LST = ['a', 'b', 'c']
# MET_REF = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
# MET_LAB = [1, 2, 3]
# MET_RAB = [1, 2, 3, 4, 5, 6, 7, 8]
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


def data_to_df(dico):
    df = pd.DataFrame(data=dico.values(), index=dico.keys())
    print(df.T)


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

    @test_for(get_fig_parameters)
    def test_get_fig_parameters(self):
        data = get_fig_parameters(classes_abondance=MET_RAB_D, parent_dict=MC_ONTO,
                                  children_dict=MC_ONTO_CH, root_item=ROOTS[METACYC],
                                  subset_abundance=MET_LAB_D, full=True, names=None)
        data_to_df(data)
