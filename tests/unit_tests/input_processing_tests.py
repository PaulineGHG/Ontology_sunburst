import unittest
from unittest.mock import patch
import io
import sys
from functools import wraps
from ontosunburst.ontosunburst import *

"""
Tests manually good file creation.
No automatic tests integrated.
"""

# ==================================================================================================
# GLOBAL
# ==================================================================================================

# --------------------------------------------------------------------------------------------------

# GLOBAL VALUES
G_LST = ['a', 'b', 'c']
G_REF = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
G_LAB = [23, 20, 5]
G_RAB = [14, 26, 20, 10, 20, 5, 4, 3, 1]
G_ONTO = {'a': ['i'], 'b': ['i'], 'c': ['j', 'k'], 'd': ['j'], 'e': ['j', 'l'],
          'f': ['k'], 'g': ['m', 'l'], 'h': ['m'],
          'i': ['x'], 'j': ['n', 'o'], 'k': ['n'],
          'l': ['x', 'o'], 'm': ['x'],
          'n': ['x'], 'o': ['v'], 'v': ['w'], 'w': ['x'], 'x': ['r'],
          'p': ['i'], 'q': ['p'], 's': ['q'], 't': ['p'], 'u': ['t']}
G_LABELS = {'r': 'Root', 'v': 'V', 'o': 'O', 'n': 'N', 'm': 'M',
            'l': 'L', 'j': 'J', 'k': 'K', 'h': 'H', 'g': 'G', 'f': 'F', 'e': 'E', 'd': 'D',
            'c': 'C', 'i': 'I', 'b': 'B'}


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

# TESTS REDUCE DAG FUNCTIONS
# --------------------------------------------------------------------------------------------------
class TestReduceDAG(unittest.TestCase):

    @test_for(get_id_to_label_dict)
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_get_id_to_label_dict(self, mock_stdout):
        id_to_label = get_id_to_label_dict(id_to_label_input=G_LABELS)
        output = mock_stdout.getvalue().strip()
        self.assertEqual(output, '3 concepts to classify\n'
                                 '3/3 concepts classified')
        self.assertEqual(classified_concepts, {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf']})