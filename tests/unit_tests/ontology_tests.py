import unittest
from unittest.mock import patch
import io
import sys
from functools import wraps
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
# EC_LST = ['1.4.5.6', '1.4.6.7', '2.1.2.3', '1.5.3', '1.6.9.-', '1.4.-.-']
EC_ONTO = {'1.4.5.-': ['1.4.-.-'], '1.4.6.-': ['1.4.-.-'], '2.1.2.-': ['2.1.-.-'],
           '1.5.3.-': ['1.5.-.-'], '1.6.9.-': ['1.6.-.-'],
           '1.4.-.-': ['1.-.-.-'], '2.1.-.-': ['2.-.-.-'], '1.5.-.-': ['1.-.-.-'],
           '1.6.-.-': ['1.-.-.-'], '1.-.-.-': [ROOTS[EC]], '2.-.-.-': [ROOTS[EC]]}
EC_ONTO_FULL = {'1.4.5.-': ['1.4.-.-'], '1.4.6.-': ['1.4.-.-'], '2.1.2.-': ['2.1.-.-'],
                '1.5.3.-': ['1.5.-.-'], '1.6.9.-': ['1.6.-.-'], '1.4.-.-': ['1.-.-.-'],
                '2.1.-.-': ['2.-.-.-'], '1.5.-.-': ['1.-.-.-'], '1.6.-.-': ['1.-.-.-'],
                '1.-.-.-': ['Enzyme'], '2.-.-.-': ['Enzyme'], '1.4.5.6': ['1.4.5.-'],
                '1.4.6.7': ['1.4.6.-'], '2.1.2.3': ['2.1.2.-'], '1.5.3': ['1.5.-.-']}

# CHEBI
# --------------------------------------------------------------------------------------------------
CH_URL = 'http://localhost:3030/chebi/'
CH_LST = ['38028', '28604', '85146']
REF_CH = ['38028', '28604', '85146',
          '23066', '27803', '37565',
          '58215', '79983', '42639']

# GO
# --------------------------------------------------------------------------------------------------
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

# TEST CLASSES EXTRACTION
# --------------------------------------------------------------------------------------------------
class TestClassesExtraction(unittest.TestCase):
    @test_for(extract_met_classes)
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_met_classes_input_ok(self, mock_stdout):
        d_obj = extract_met_classes(MET_LST, MC_ONTO)
        output = mock_stdout.getvalue().strip()
        self.assertEqual(output, '3 metabolic objects to classify\n'
                                 '3/3 metabolic objects classified')
        self.assertEqual(d_obj, {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf']})

    @test_for(extract_met_classes)
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_met_classes_input_errors(self, mock_stdout):
        d_obj = extract_met_classes(MET_LST + ['x'], MC_ONTO)
        output = mock_stdout.getvalue().strip()
        self.assertEqual(output, '4 metabolic objects to classify\n'
                                 'x not classified.\n'
                                 '3/4 metabolic objects classified')
        self.assertEqual(d_obj, {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf']})

    @test_for(extract_ec_classes)
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_ec_classes_input_ok(self, mock_stdout):
        d_obj, d_onto = extract_ec_classes(EC_LST, EC_ONTO)
        output = mock_stdout.getvalue().strip()
        wanted_d_obj = {'1.4.5.6': ['1.4.5.-'], '1.4.6.7': ['1.4.6.-'], '2.1.2.3': ['2.1.2.-'],
                        '1.5.3': ['1.5.-.-'], '1.6.9.-': ['1.6.-.-'], '1.-.-.-': ['Enzyme'],
                        '1.4.-.-': ['1.-.-.-']}
        self.assertEqual(output, '7 EC numbers to classify\n'
                                 '7/7 EC numbers classified')
        self.assertEqual(d_obj, wanted_d_obj)
        self.assertEqual(d_onto, EC_ONTO_FULL)

    @test_for(extract_ec_classes)
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_ec_classes_input_errors(self, mock_stdout):
        d_obj, d_onto = extract_ec_classes(EC_LST + ['3.5.6.9', 'ecID'], EC_ONTO)
        output = mock_stdout.getvalue().strip()
        wanted_d_obj = {'1.4.5.6': ['1.4.5.-'], '1.4.6.7': ['1.4.6.-'], '2.1.2.3': ['2.1.2.-'],
                        '1.5.3': ['1.5.-.-'], '1.6.9.-': ['1.6.-.-'], '1.-.-.-': ['Enzyme'],
                        '1.4.-.-': ['1.-.-.-']}
        self.assertEqual(output, '9 EC numbers to classify\n'
                                 '3.5.6.9 not classified\n'
                                 'ecID not classified\n'
                                 '7/9 EC numbers classified')
        self.assertEqual(d_obj, wanted_d_obj)
        self.assertEqual(d_onto, EC_ONTO_FULL)

    @test_for(get_parents)
    def test_get_parents_linear_path(self):
        # Simple linear direction
        parents = get_parents('a', {'ab'}, MC_ONTO, ROOTS[METACYC])
        self.assertEqual(parents, {'FRAMES', 'ab'})

    @test_for(get_parents)
    def test_get_parents_complex_path(self):
        # With multiple parents having multiple parents and different size of path until root
        parents = get_parents('c', {'cde', 'cf'}, MC_ONTO, ROOTS[METACYC])
        self.assertEqual(parents, {'cdeeg+', 'FRAMES', 'cf', 'cde', 'cdecf', 'cdeeg'})

    @test_for(get_parents)
    def test_get_parents_ec(self):
        # With EC (simple path)
        parents = get_parents('2.1.2.3', {'2.1.2.-'}, EC_ONTO_FULL, ROOTS[EC])
        self.assertEqual(parents, {'Enzyme', '2.1.-.-', '2.-.-.-', '2.1.2.-'})

    @test_for(get_parents)
    def test_get_parents_direct_root_child(self):
        # With direct child of root item
        parents = get_parents('1.-.-.-', {'Enzyme'}, EC_ONTO_FULL, ROOTS[EC])
        self.assertEqual(parents, {'Enzyme'})

    @test_for(get_all_classes)
    def test_get_all_classes(self):
        leaf_classes = {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf']}
        all_classes_met = get_all_classes(leaf_classes, MC_ONTO, ROOTS[METACYC])
        wanted_all_classes = {'a': {'FRAMES', 'ab'}, 'b': {'FRAMES', 'ab'},
                              'c': {'cdeeg+', 'cde', 'cdeeg', 'FRAMES', 'cdecf', 'cf'}}
        self.assertEqual(all_classes_met, wanted_all_classes)

    @test_for(get_all_classes)
    def test_get_all_classes_ec(self):
        ec_leaf_classes = {'1.4.5.6': ['1.4.5.-'], '1.4.6.7': ['1.4.6.-'], '2.1.2.3': ['2.1.2.-'],
                           '1.5.3': ['1.5.-.-'], '1.6.9.-': ['1.6.-.-'], '1.-.-.-': ['Enzyme'],
                           '1.4.-.-': ['1.-.-.-']}
        wanted_all_classes = {'1.4.5.6': {'1.4.5.-', '1.4.-.-', 'Enzyme', '1.-.-.-'},
                              '1.4.6.7': {'Enzyme', '1.4.-.-', '1.4.6.-', '1.-.-.-'},
                              '2.1.2.3': {'2.-.-.-', '2.1.-.-', '2.1.2.-', 'Enzyme'},
                              '1.5.3': {'1.5.-.-', 'Enzyme', '1.-.-.-'},
                              '1.6.9.-': {'1.6.-.-', 'Enzyme', '1.-.-.-'},
                              '1.-.-.-': {'Enzyme'}, '1.4.-.-': {'Enzyme', '1.-.-.-'}}
        all_classes_ec = get_all_classes(ec_leaf_classes, EC_ONTO, ROOTS[EC])
        self.assertEqual(all_classes_ec, wanted_all_classes)

    @test_for(extract_classes)
    def test_extract_classes_metacyc(self):
        mc_classes, d_classes_ontology = extract_classes(METACYC, MET_LST, ROOTS[METACYC], MC_ONTO,
                                                         None)
        wanted_mc_classes = {'a': {'FRAMES', 'ab'}, 'b': {'FRAMES', 'ab'},
                             'c': {'cdeeg+', 'cde', 'cdeeg', 'FRAMES', 'cdecf', 'cf'}}
        self.assertEqual(mc_classes, wanted_mc_classes)
        self.assertTrue(dicts_with_sorted_lists_equal(d_classes_ontology, MC_ONTO))

    @test_for(extract_classes)
    def test_extract_classes_kegg(self):
        kg_classes, d_classes_ontology = extract_classes(KEGG, MET_LST, ROOTS[KEGG], KG_ONTO, None)
        print(kg_classes)
        wanted_kg_classes = {'a': {'ab', 'kegg'}, 'b': {'ab', 'kegg'},
                             'c': {'cde', 'cdeeg+', 'kegg', 'cdecf', 'cdeeg', 'cf'}}
        self.assertEqual(kg_classes, wanted_kg_classes)
        self.assertTrue(dicts_with_sorted_lists_equal(d_classes_ontology, KG_ONTO))

    @test_for(extract_classes)
    def test_extract_classes_ec(self):
        ec_classes, d_classes_ontology = extract_classes(EC, EC_LST, ROOTS[EC], EC_ONTO, None)
        print(ec_classes)
        wanted_ec_classes = {'1.4.5.6': {'1.4.5.-', '1.4.-.-', 'Enzyme', '1.-.-.-'},
                             '1.4.6.7': {'Enzyme', '1.4.-.-', '1.4.6.-', '1.-.-.-'},
                             '2.1.2.3': {'2.-.-.-', '2.1.-.-', '2.1.2.-', 'Enzyme'},
                             '1.5.3': {'1.5.-.-', 'Enzyme', '1.-.-.-'},
                             '1.6.9.-': {'1.6.-.-', 'Enzyme', '1.-.-.-'},
                             '1.-.-.-': {'Enzyme'}, '1.4.-.-': {'Enzyme', '1.-.-.-'}}
        self.assertEqual(ec_classes, wanted_ec_classes)
        self.assertTrue(dicts_with_sorted_lists_equal(d_classes_ontology, EC_ONTO_FULL))


# TEST CHEBI + GO : NEED DEDICATED SPARQL SERVER RUNNING TO EXECUTE
# --------------------------------------------------------------------------------------------------
# ChEBI
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


# GO
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


# TEST ABUNDANCES
# --------------------------------------------------------------------------------------------------

class TestAbundances(unittest.TestCase):
    @test_for(get_abundance_dict)
    def test_get_abundance_dict_abundances_not_ref(self):
        abundance_dict = get_abundance_dict(abundances=MET_LAB, metabolic_objects=MET_LST,
                                            ref=False)
        self.assertEqual(abundance_dict, {'a': 1, 'b': 2, 'c': 3})

    @test_for(get_abundance_dict)
    def test_get_abundance_dict_no_abundances_not_ref(self):
        abundance_dict = get_abundance_dict(abundances=None, metabolic_objects=MET_LST,
                                            ref=False)
        self.assertEqual(abundance_dict, {'a': 1, 'b': 1, 'c': 1})

    @test_for(get_abundance_dict)
    def test_get_abundance_dict_abundances_ref(self):
        abundance_dict = get_abundance_dict(abundances=MET_RAB, metabolic_objects=MET_REF,
                                            ref=True)
        self.assertEqual(abundance_dict, {'a': 1, 'b': 2, 'c': 3, 'd': 4,
                                          'e': 5, 'f': 6, 'g': 7, 'h': 8})

    @test_for(get_abundance_dict)
    def test_get_abundance_dict_no_abundances_ref(self):
        abundance_dict = get_abundance_dict(abundances=None, metabolic_objects=MET_REF,
                                            ref=True)
        self.assertEqual(abundance_dict, {'a': 1, 'b': 1, 'c': 1, 'd': 1,
                                          'e': 1, 'f': 1, 'g': 1, 'h': 1})

    @test_for(get_abundance_dict)
    def test_get_abundance_dict_errors_not_ref(self):
        with self.assertRaises(AttributeError) as e:
            get_abundance_dict(abundances=MET_LAB + [4], metabolic_objects=MET_LST, ref=False)
        self.assertEqual(str(e.exception), 'Length of "metabolic_objects" parameter must be '
                                           'equal to "abundances" parameter length : 3 != 4')

    @test_for(get_abundance_dict)
    def test_get_abundance_dict_errors_ref(self):
        with self.assertRaises(AttributeError) as e:
            get_abundance_dict(abundances=MET_RAB[:-1], metabolic_objects=MET_REF, ref=True)
        self.assertEqual(str(e.exception), 'Length of "reference_set" parameter must be '
                                           'equal to "ref_abundances" parameter length : 8 != 7')

    @test_for(get_classes_abundance)
    def test_get_classes_abundance(self):
        classes_abundances = get_classes_abundance(all_classes={}, abundances_dict={},
                                                   show_leaves=False)
