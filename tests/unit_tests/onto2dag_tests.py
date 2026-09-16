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
        sub_dag = SubDAG(I_DCT, R_DCT, ALL_CPT, ONTO, ROOT, LABELS)
        for node in sub_dag.nodes:
            print(node.onto_id, node.label, node.experimental_weight, node.cumulative_weight,
                  node.proportion, node.ref_experimental_weight, node.ref_cumulative_weight,
                  node.ref_proportion, node.parents, node.children, node.difference,
                  node.enrichment_p_val)

# TESTS REDUCE DAG FUNCTIONS
# --------------------------------------------------------------------------------------------------
# class TestReduceDAG(unittest.TestCase):
#
#     @test_for(classify_concepts)
#     @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
#     def test_classify_concepts_ok(self, mock_stdout):
#         classified_concepts = classify_concepts(concepts=CPT_LST, ontology_dag=ONTO_DAG)
#         output = mock_stdout.getvalue().strip()
#         self.assertEqual(output, '3 concepts to classify\n'
#                                  '3/3 concepts classified')
#         self.assertEqual(classified_concepts, {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf']})
#
#     @test_for(classify_concepts)
#     @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
#     def test_classify_concepts_errors(self, mock_stdout):
#         classified_concepts = classify_concepts(concepts=CPT_LST + ['x'], ontology_dag=ONTO_DAG)
#         output = mock_stdout.getvalue().strip()
#         self.assertEqual(output, '4 concepts to classify\n'
#                                  'x not classified.\n'
#                                  '3/4 concepts classified')
#         self.assertEqual(classified_concepts, {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf']})
#
#     @test_for(get_parents)
#     def test_get_parents_linear_path(self):
#         # Simple linear direction
#         parents = get_parents('a', {'ab'}, ONTO_DAG, ROOT)
#         self.assertEqual(parents, {'root', 'ab'})
#
#     @test_for(get_parents)
#     def test_get_parents_complex_path(self):
#         # With multiple parents having multiple parents and different size of path until root
#         parents = get_parents('c', {'cde', 'cf'}, ONTO_DAG, ROOT)
#         self.assertEqual(parents, {'cdeeg+', 'root', 'cf', 'cde', 'cdecf', 'cdeeg'})
#
#     @test_for(get_all_classes)
#     def test_get_all_classes(self):
#         leaf_classes = {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf']}
#         all_classes_met = get_all_classes(leaf_classes, ONTO_DAG, ROOT)
#         wanted_all_classes = {'a': {'root', 'ab'}, 'b': {'root', 'ab'},
#                               'c': {'cdeeg+', 'cde', 'cdeeg', 'root', 'cdecf', 'cf'}}
#         self.assertEqual(all_classes_met, wanted_all_classes)
#
#
# # TESTS WEIGHTS CALCULATION
# # --------------------------------------------------------------------------------------------------
# class TestWeightsCalculation(unittest.TestCase):
#     @test_for(get_abundance_dict)
#     def test_get_abundance_dict_abundances(self):
#         abundance_dict = get_abundance_dict(abundances=CPT_AB, concepts=CPT_LST)
#         self.assertEqual(abundance_dict, {'a': 1, 'b': 2, 'c': 3})
#
#     @test_for(get_abundance_dict)
#     def test_get_abundance_dict_no_abundances(self):
#         abundance_dict = get_abundance_dict(abundances=None, concepts=CPT_LST)
#         self.assertEqual(abundance_dict, {'a': 1, 'b': 1, 'c': 1})
#
#     @test_for(get_abundance_dict)
#     def test_get_abundance_dict_abundances_ref(self):
#         abundance_dict = get_abundance_dict(abundances=RCPT_AB, concepts=RCPT_LST)
#         self.assertEqual(abundance_dict, {'a': 1, 'b': 2, 'c': 3, 'd': 4,
#                                           'e': 5, 'f': 6, 'g': 7, 'h': 8})
#
#     @test_for(get_abundance_dict)
#     def test_get_abundance_dict_no_abundances_ref(self):
#         abundance_dict = get_abundance_dict(abundances=None, concepts=RCPT_LST)
#         self.assertEqual(abundance_dict, {'a': 1, 'b': 1, 'c': 1, 'd': 1,
#                                           'e': 1, 'f': 1, 'g': 1, 'h': 1})
#
#     @test_for(get_abundance_dict)
#     def test_get_abundance_dict_errors(self):
#         with self.assertRaises(AttributeError) as e:
#             get_abundance_dict(abundances=CPT_AB + [4], concepts=CPT_LST)
#         self.assertEqual(str(e.exception), 'Length of concepts IDs list must '
#                                            'be equal to its abundances list length : 3 != 4')
#
#     @test_for(calculate_weights)
#     def test_get_classes_abundance_leaves(self):
#         all_classes = {'a': {'root', 'ab'}, 'b': {'root', 'ab'},
#                        'c': {'cdecf', 'cdeeg+', 'root', 'cde', 'cdeeg', 'cf'},
#                        'd': {'cdecf', 'cdeeg+', 'root', 'cde', 'cdeeg'},
#                        'e': {'cdeeg+', 'root', 'cde', 'cdecf', 'eg', 'cdeeg'},
#                        'f': {'cdecf', 'root', 'cf'},
#                        'g': {'cdeeg', 'cdeeg+', 'root', 'eg', 'gh'}, 'h': {'root', 'gh'}}
#         abundances_dict = {'a': 1, 'b': 2, 'c': 3, 'd': 4, 'e': 5, 'f': 6, 'g': 7, 'h': 8}
#         classes_abundances = calculate_weights(all_classes, abundances_dict, show_leaves=True)
#         wanted_abundances = {'root': 36, 'cdeeg+': 19, 'cdeeg': 19, 'cdecf': 18, 'gh': 15,
#                              'eg': 12, 'cde': 12, 'cf': 9, 'h': 8, 'g': 7, 'f': 6, 'e': 5,
#                              'd': 4, 'c': 3, 'ab': 3, 'b': 2, 'a': 1}
#         self.assertEqual(classes_abundances, wanted_abundances)
#
#     @test_for(calculate_weights)
#     def test_get_classes_abundance_leaves_sub(self):
#         all_classes = {'a': {'root', 'ab'}, 'b': {'root', 'ab'},
#                        'c': {'cdeeg+', 'cde', 'cdeeg', 'root', 'cdecf', 'cf'}}
#         abundances_dict = {'a': 1, 'b': 2, 'c': 3}
#         classes_abundances = calculate_weights(all_classes, abundances_dict, show_leaves=True)
#         wanted_abundances = {'root': 6, 'cde': 3, 'cf': 3, 'cdecf': 3, 'cdeeg+': 3, 'cdeeg': 3,
#                              'c': 3, 'ab': 3, 'b': 2, 'a': 1}
#         self.assertEqual(classes_abundances, wanted_abundances)
#
#     @test_for(calculate_weights)
#     def test_get_classes_abundance_no_leaves(self):
#         all_classes = {'a': {'root', 'ab'}, 'b': {'root', 'ab'},
#                        'c': {'cdecf', 'cdeeg+', 'root', 'cde', 'cdeeg', 'cf'},
#                        'd': {'cdecf', 'cdeeg+', 'root', 'cde', 'cdeeg'},
#                        'e': {'cdeeg+', 'root', 'cde', 'cdecf', 'eg', 'cdeeg'},
#                        'f': {'cdecf', 'root', 'cf'},
#                        'g': {'cdeeg', 'cdeeg+', 'root', 'eg', 'gh'}, 'h': {'root', 'gh'}}
#         abundances_dict = {'a': 1, 'b': 2, 'c': 3, 'd': 4, 'e': 5, 'f': 6, 'g': 7, 'h': 8}
#         classes_abundances = calculate_weights(all_classes, abundances_dict, show_leaves=False)
#         wanted_abundances = {'root': 36, 'cdeeg+': 19, 'cdeeg': 19, 'cdecf': 18, 'gh': 15,
#                              'eg': 12, 'cde': 12, 'cf': 9, 'ab': 3}
#         self.assertEqual(classes_abundances, wanted_abundances)
#
#     @test_for(calculate_weights)
#     def test_get_classes_abundance_no_leaves_sub(self):
#         all_classes = {'a': {'root', 'ab'}, 'b': {'root', 'ab'},
#                        'c': {'cdeeg+', 'cde', 'cdeeg', 'root', 'cdecf', 'cf'}}
#         abundances_dict = {'a': 1, 'b': 2, 'c': 3}
#         classes_abundances = calculate_weights(all_classes, abundances_dict, show_leaves=False)
#         wanted_abundances = {'root': 6, 'cde': 3, 'cf': 3, 'cdecf': 3, 'cdeeg+': 3, 'cdeeg': 3,
#                              'ab': 3}
#         self.assertEqual(classes_abundances, wanted_abundances)
#
#     @test_for(calculate_weights)
#     def test_get_classes_abundance_different_level_abundances(self):
#         all_classes = {'c': {'cdecf', 'cdeeg+', 'root', 'cde', 'cdeeg', 'cf'},
#                        'd': {'cdecf', 'cdeeg+', 'root', 'cde', 'cdeeg'},
#                        'e': {'cdeeg+', 'root', 'cde', 'cdecf', 'eg', 'cdeeg'},
#                        'f': {'cdecf', 'root', 'cf'}, 'cf': {'cdecf', 'root'}}
#         abundances_dict = {'c': 3, 'd': 4, 'e': 5, 'f': 2, 'cf': 2}
#         classes_abundances = calculate_weights(all_classes, abundances_dict, show_leaves=True)
#         wanted_abundances = {ROOT: 16, 'cdecf': 16, 'cde': 12, 'cdeeg+': 12, 'cdeeg': 12,
#                              'cf': 7, 'eg': 5, 'e': 5, 'd': 4, 'c': 3, 'f': 2}
#         self.assertEqual(classes_abundances, wanted_abundances)
#
#
# # TEST MAIN FUNCTIONS
# # --------------------------------------------------------------------------------------------------
#
# class TestMainFunctions(unittest.TestCase):
#
#     @test_for(reduce_d_ontology)
#     def test_reduce_d_ontology_onto_dag(self):
#         classes_abundance = {'root': 6, 'cde': 3, 'cf': 3, 'cdecf': 3, 'cdeeg+': 3, 'cdeeg': 3,
#                              'ab': 3}
#         d_ontology_reduced = reduce_d_ontology(ONTO_DAG, classes_abundance)
#         wanted_d_ontology_reduced = {'ab': ['root'], 'cde': ['cdecf', 'cdeeg'],
#                                      'cf': ['cdecf'], 'cdecf': ['root'], 'cdeeg': ['cdeeg+'],
#                                      'cdeeg+': ['root']}
#         self.assertEqual(d_ontology_reduced, wanted_d_ontology_reduced)
#
#     @test_for(reduce_d_ontology)
#     def test_reduce_d_ontology_id_to_label(self):
#         classes_abundance = {'root': 6, 'cde': 3, 'cf': 3, 'cdecf': 3, 'cdeeg+': 3, 'cdeeg': 3,
#                              'ab': 3}
#         d_ontology_reduced = reduce_d_ontology(ID2LAB, classes_abundance)
#         wanted_d_ontology_reduced = {'root': 'Root', 'cdeeg+': 'CDEEG+', 'cdeeg': 'CDEEG',
#                                      'cdecf': 'CDECF', 'cde': 'CDE', 'cf': 'CF', 'ab': 'AB'}
#         self.assertEqual(d_ontology_reduced, wanted_d_ontology_reduced)
#
#     @test_for(ontology_to_weighted_dag)
#     def test_ontology_to_weighted_dag_no_lvs(self):
#         calculated_weights = ontology_to_weighted_dag(CPT_LST, CPT_AB, ROOT, ONTO_DAG, False)
#         wanted_abundances = {'root': 6, 'cde': 3, 'cf': 3, 'cdecf': 3, 'cdeeg+': 3, 'cdeeg': 3,
#                              'ab': 3}
#         self.assertEqual(calculated_weights, wanted_abundances)
#
#     @test_for(ontology_to_weighted_dag)
#     def test_ontology_to_weighted_dag_ref_lvs(self):
#         calculated_weights = ontology_to_weighted_dag(RCPT_LST, RCPT_AB, ROOT, ONTO_DAG, True)
#         wanted_abundances = {'root': 36, 'cdeeg+': 19, 'cdeeg': 19, 'cdecf': 18, 'gh': 15,
#                              'eg': 12, 'cde': 12, 'cf': 9, 'h': 8, 'g': 7, 'f': 6, 'e': 5,
#                              'd': 4, 'c': 3, 'ab': 3, 'b': 2, 'a': 1}
#         self.assertEqual(calculated_weights, wanted_abundances)
