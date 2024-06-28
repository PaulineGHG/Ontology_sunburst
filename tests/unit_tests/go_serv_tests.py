import unittest
from unittest.mock import patch
import io
import sys
from functools import wraps
from ontosunburst.ontology import *

# ==================================================================================================
# GLOBAL
# ==================================================================================================

GO_LST = ['GO:0043227', 'GO:0043229', 'GO:0043231', 'GO:0044422']
GO_URL = 'http://localhost:3030/go/'


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

# TEST ONTOLOGY : CLASSES EXTRACTION
# --------------------------------------------------------------------------------------------------
class TestGOClassesExtraction(unittest.TestCase):
    @test_for(extract_go_classes)
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_go_classes(self, mock_stdout):
        all_classes, d_classes_ontology = extract_go_classes(GO_LST, GO_URL)
        output = mock_stdout.getvalue().strip()
        wanted_ontology = {'membrane-bounded organelle': ['organelle'],
                           'organelle': ['cellular anatomical entity'],
                           'cellular_component': ['GO'],
                           'cellular anatomical entity': ['cellular_component'],
                           'go:0043227': ['membrane-bounded organelle'],
                           'intracellular organelle': ['organelle'],
                           'go:0043229': ['intracellular organelle'],
                           'intracellular membrane-bounded organelle':
                               ['membrane-bounded organelle', 'intracellular organelle'],
                           'go:0043231': ['intracellular membrane-bounded organelle']}
        wanted_classes = {'go:0043227': {'organelle', 'cellular anatomical entity',
                                         'membrane-bounded organelle', 'cellular_component', 'GO'},
                          'go:0043229': {'organelle', 'cellular anatomical entity',
                                         'intracellular organelle', 'cellular_component', 'GO'},
                          'go:0043231': {'organelle', 'cellular anatomical entity',
                                         'intracellular organelle', 'membrane-bounded organelle',
                                         'cellular_component',
                                         'intracellular membrane-bounded organelle', 'GO'}}

        self.assertEqual(output, 'No GO class found for : go:0044422')
        self.assertDictEqual(all_classes, wanted_classes)
        self.assertTrue(dicts_with_sorted_lists_equal(d_classes_ontology, wanted_ontology))

    @test_for(extract_classes)
    def test_extract_classes_go(self):
        go_classes, d_classes_ontology = extract_classes(GO, GO_LST, ROOTS[GO], None, GO_URL)
        wanted_ontology = {'membrane-bounded organelle': ['organelle'],
                           'organelle': ['cellular anatomical entity'],
                           'cellular_component': ['GO'],
                           'cellular anatomical entity': ['cellular_component'],
                           'go:0043227': ['membrane-bounded organelle'],
                           'intracellular organelle': ['organelle'],
                           'go:0043229': ['intracellular organelle'],
                           'intracellular membrane-bounded organelle':
                               ['membrane-bounded organelle', 'intracellular organelle'],
                           'go:0043231': ['intracellular membrane-bounded organelle']}
        wanted_classes = {'go:0043227': {'organelle', 'cellular anatomical entity',
                                         'membrane-bounded organelle', 'cellular_component', 'GO'},
                          'go:0043229': {'organelle', 'cellular anatomical entity',
                                         'intracellular organelle', 'cellular_component', 'GO'},
                          'go:0043231': {'organelle', 'cellular anatomical entity',
                                         'intracellular organelle', 'membrane-bounded organelle',
                                         'cellular_component',
                                         'intracellular membrane-bounded organelle', 'GO'}}
        self.assertEqual(go_classes, wanted_classes)
        self.assertTrue(dicts_with_sorted_lists_equal(d_classes_ontology, wanted_ontology))
