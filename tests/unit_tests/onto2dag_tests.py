import unittest
from unittest.mock import patch
import io
import sys
from functools import wraps
from ontosunburst.onto2dag import *

"""
Tests manually good file creation.
No automatic tests integrated.
"""

# ==================================================================================================
# GLOBAL
# ==================================================================================================

# --------------------------------------------------------------------------------------------------

# GLOBAL VALUES
I_LST = ['a', 'b', 'c']
R_LST = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
I_DCT = {'a': 23, 'b': 20, 'c': 5}
R_DCT = {'a': 14, 'b': 26, 'c': 20, 'd': 10, 'e': 20, 'f': 5, 'g': 4, 'h': 3, 'i': 1}
ALL_CPT = set(R_LST)
ONTO = {'a': ['i'], 'b': ['i'], 'c': ['j', 'k'], 'd': ['j'], 'e': ['j', 'l'],
        'f': ['k'], 'g': ['m', 'l'], 'h': ['m'],
        'i': ['x'], 'j': ['n', 'o'], 'k': ['n'],
        'l': ['x', 'o'], 'm': ['x'],
        'n': ['x'], 'o': ['v'], 'v': ['w'], 'w': ['x'], 'x': ['r'],
        'p': ['i'], 'q': ['p'], 's': ['q'], 't': ['p'], 'u': ['t']}
LABELS = {'r': 'Root', 'v': 'V', 'o': 'O', 'n': 'N', 'm': 'M',
          'l': 'L', 'j': 'J', 'k': 'K', 'h': 'H', 'g': 'G', 'f': 'F', 'e': 'E', 'd': 'D',
          'c': 'C', 'i': 'I', 'b': 'B'}
ROOT = 'r'
ANCESTORS = {'a': {'x', 'r', 'i'},
             'b': {'x', 'r', 'i'},
             'c': {'v', 'r', 'w', 'o', 'j', 'x', 'k', 'n'},
             'd': {'v', 'r', 'w', 'o', 'j', 'x', 'n'},
             'e': {'v', 'r', 'w', 'o', 'j', 'x', 'n', 'l'},
             'f': {'n', 'x', 'k', 'r'},
             'g': {'v', 'r', 'w', 'm', 'o', 'x', 'l'},
             'h': {'x', 'm', 'r'},
             'i': {'x', 'r'}}
CML_W = {'x': 48, 'r': 48, 'i': 43, 'a': 23, 'b': 20, 'j': 5,
         'k': 5, 'v': 5, 'n': 5, 'w': 5, 'o': 5, 'c': 5}
R_CML_W = {'r': 103, 'x': 103, 'n': 55, 'w': 54, 'o': 54, 'v': 54, 'j': 50, 'i': 41, 'b': 26,
           'k': 25, 'l': 24, 'e': 20, 'c': 20, 'a': 14, 'd': 10, 'm': 7, 'f': 5, 'g': 4, 'h': 3}


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

class TestOntoToDag(unittest.TestCase):

    @test_for(get_ancestors)
    @test_for(get_ancestors_recursively)
    def test_get_ancestors(self):
        all_parents = get_ancestors(ALL_CPT, ONTO, ROOT)
        self.assertDictEqual(all_parents, ANCESTORS)

    @test_for(get_cumulative_w)
    def test_get_cumulative_w(self):
        i_cml_w = get_cumulative_w(ANCESTORS, I_DCT)
        self.assertDictEqual(i_cml_w, CML_W)

    @test_for(get_cumulative_w)
    def test_get_cumulative_w(self):
        r_cml_w = get_cumulative_w(ANCESTORS, R_DCT)
        self.assertDictEqual(r_cml_w, R_CML_W)

    def test_sub_dag(self):
        sub_dag = SubDAG(I_DCT, R_DCT, ALL_CPT, ONTO, ROOT, LABELS, BINOMIAL_TEST)
        for node in sub_dag.nodes:
            # node._print_arguments()
            print(node.onto_id, node.difference)

