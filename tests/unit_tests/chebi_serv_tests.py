import unittest
from unittest.mock import patch
import io
import sys
from functools import wraps
from ontosunburst.ontology import *

# ==================================================================================================
# GLOBAL
# ==================================================================================================

CH_URL = 'http://localhost:3030/chebi/'
CH_LST = ['38028', '28604', '85146']
REF_CH = ['38028', '28604', '85146',
          '23066', '27803', '37565',
          '58215', '79983', '42639']


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
class TestChEBIClassesExtraction(unittest.TestCase):
    @test_for(extract_chebi_roles)
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_chebi_roles(self, mock_stdout):
        all_roles, d_roles_ontology = extract_chebi_roles(CH_LST, CH_URL)
        output = mock_stdout.getvalue().strip()
        wanted_ontology = {'xenobiotic': ['biological role'], 'biological role': ['role'],
                           '38028': ['xenobiotic'], 'algal metabolite': ['eukaryotic metabolite'],
                           'eukaryotic metabolite': ['metabolite'],
                           'metabolite': ['biochemical role'],
                           'biochemical role': ['biological role'],
                           'marine metabolite': ['metabolite'],
                           'plant metabolite': ['eukaryotic metabolite'],
                           'animal metabolite': ['eukaryotic metabolite'],
                           '28604': ['animal metabolite', 'plant metabolite', 'algal metabolite',
                                     'marine metabolite'],
                           'food emulsifier': ['food additive', 'emulsifier'],
                           'food additive': ['food component', 'application'],
                           'application': ['role'], 'food component': ['physiological role'],
                           'physiological role': ['biological role'],
                           'emulsifier': ['chemical role'], 'chemical role': ['role'],
                           'laxative': ['drug'], 'drug': ['pharmaceutical'],
                           'pharmaceutical': ['application'],
                           '85146': ['food emulsifier', 'laxative']}
        wanted_roles = {'38028': {'role', 'biological role', 'xenobiotic'},
                        '28604': {'algal metabolite', 'biochemical role', 'biological role',
                                  'eukaryotic metabolite', 'marine metabolite', 'plant metabolite',
                                  'metabolite', 'animal metabolite', 'role'},
                        '85146': {'emulsifier', 'food additive', 'pharmaceutical', 'application',
                                  'biological role', 'food component', 'chemical role', 'role',
                                  'physiological role', 'food emulsifier', 'drug', 'laxative'}}

        self.assertEqual(output, '3/3 chebi id with roles associated.')
        self.assertDictEqual(all_roles, wanted_roles)
        self.assertTrue(dicts_with_sorted_lists_equal(d_roles_ontology, wanted_ontology))

    @test_for(extract_classes)
    def test_extract_classes_chebi(self):
        ch_classes, d_classes_ontology = extract_classes(CHEBI, CH_LST, ROOTS[CHEBI], None, CH_URL)
        wanted_ontology = {'xenobiotic': ['biological role'], 'biological role': ['role'],
                           '38028': ['xenobiotic'], 'algal metabolite': ['eukaryotic metabolite'],
                           'eukaryotic metabolite': ['metabolite'],
                           'metabolite': ['biochemical role'],
                           'biochemical role': ['biological role'],
                           'marine metabolite': ['metabolite'],
                           'plant metabolite': ['eukaryotic metabolite'],
                           'animal metabolite': ['eukaryotic metabolite'],
                           '28604': ['animal metabolite', 'plant metabolite', 'algal metabolite',
                                     'marine metabolite'],
                           'food emulsifier': ['food additive', 'emulsifier'],
                           'food additive': ['food component', 'application'],
                           'application': ['role'], 'food component': ['physiological role'],
                           'physiological role': ['biological role'],
                           'emulsifier': ['chemical role'], 'chemical role': ['role'],
                           'laxative': ['drug'], 'drug': ['pharmaceutical'],
                           'pharmaceutical': ['application'],
                           '85146': ['food emulsifier', 'laxative']}
        wanted_roles = {'38028': {'role', 'biological role', 'xenobiotic'},
                        '28604': {'algal metabolite', 'biochemical role', 'biological role',
                                  'eukaryotic metabolite', 'marine metabolite', 'plant metabolite',
                                  'metabolite', 'animal metabolite', 'role'},
                        '85146': {'emulsifier', 'food additive', 'pharmaceutical', 'application',
                                  'biological role', 'food component', 'chemical role', 'role',
                                  'physiological role', 'food emulsifier', 'drug', 'laxative'}}
        self.assertEqual(ch_classes, wanted_roles)
        self.assertTrue(dicts_with_sorted_lists_equal(d_classes_ontology, wanted_ontology))
