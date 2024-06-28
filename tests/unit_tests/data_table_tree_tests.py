import unittest
from unittest.mock import patch
import io
import sys
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

MET_LST = ['a', 'b', 'c']
MET_REF = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
MET_LAB = [1, 2, 3]
MET_RAB = [1, 2, 3, 4, 5, 6, 7, 8]
MC_ONTO = {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf'], 'd': ['cde'], 'e': ['cde', 'eg'],
           'f': ['cf'], 'g': ['gh', 'eg'], 'h': ['gh'],
           'ab': [ROOTS[METACYC]], 'cde': ['cdecf', 'cdeeg'], 'cf': ['cdecf'],
           'eg': [ROOTS[METACYC], 'cdeeg'], 'gh': [ROOTS[METACYC]],
           'cdecf': [ROOTS[METACYC]], 'cdeeg': ['cdeeg+'], 'cdeeg+': [ROOTS[METACYC]]}
KG_ONTO = {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf'], 'd': ['cde'], 'e': ['cde', 'eg'],
           'f': ['cf'], 'g': ['gh', 'eg'], 'h': ['gh'],
           'ab': [ROOTS[KEGG]], 'cde': ['cdecf', 'cdeeg'], 'cf': ['cdecf'],
           'eg': [ROOTS[KEGG], 'cdeeg'], 'gh': [ROOTS[KEGG]],
           'cdecf': [ROOTS[KEGG]], 'cdeeg': ['cdeeg+'], 'cdeeg+': [ROOTS[KEGG]]}

# EC
# --------------------------------------------------------------------------------------------------
EC_LST = ['1.4.5.6', '1.4.6.7', '2.1.2.3', '1.5.3', '1.6.9.-', '1.-.-.-', '1.4.-.-']
EC_ONTO = {'1.4.5.-': ['1.4.-.-'], '1.4.6.-': ['1.4.-.-'], '2.1.2.-': ['2.1.-.-'],
           '1.5.3.-': ['1.5.-.-'], '1.6.9.-': ['1.6.-.-'],
           '1.4.-.-': ['1.-.-.-'], '2.1.-.-': ['2.-.-.-'], '1.5.-.-': ['1.-.-.-'],
           '1.6.-.-': ['1.-.-.-'], '1.-.-.-': [ROOTS[EC]], '2.-.-.-': [ROOTS[EC]]}
EC_ONTO_FULL = {'1.4.5.-': ['1.4.-.-'], '1.4.6.-': ['1.4.-.-'], '2.1.2.-': ['2.1.-.-'],
                '1.5.3.-': ['1.5.-.-'], '1.6.9.-': ['1.6.-.-'], '1.4.-.-': ['1.-.-.-'],
                '2.1.-.-': ['2.-.-.-'], '1.5.-.-': ['1.-.-.-'], '1.6.-.-': ['1.-.-.-'],
                '1.-.-.-': [ROOTS[EC]], '2.-.-.-': [ROOTS[EC]], '1.4.5.6': ['1.4.5.-'],
                '1.4.6.7': ['1.4.6.-'], '2.1.2.3': ['2.1.2.-'], '1.5.3': ['1.5.-.-']}


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

class TestDataTable(unittest.TestCase):
    @test_for(extract_met_classes)
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_met_classes_input_ok(self, mock_stdout):
        output = mock_stdout.getvalue().strip()
        self.assertEqual(output, '')
