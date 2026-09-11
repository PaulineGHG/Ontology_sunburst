import copy
import json
import os.path
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
I_LST = ['a', 'b', 'c']
R_LST = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
I_DCT = {'a': 23, 'b': 20, 'c': 5}
R_DCT = {'a': 14, 'b': 26, 'c': 20, 'd': 10, 'e': 20, 'f': 5, 'g': 4, 'h': 3, 'i': 1}
ONTO = {'a': ['i'], 'b': ['i'], 'c': ['j', 'k'], 'd': ['j'], 'e': ['j', 'l'],
        'f': ['k'], 'g': ['m', 'l'], 'h': ['m'],
        'i': ['x'], 'j': ['n', 'o'], 'k': ['n'],
        'l': ['x', 'o'], 'm': ['x'],
        'n': ['x'], 'o': ['v'], 'v': ['w'], 'w': ['x'], 'x': ['r'],
        'p': ['i'], 'q': ['p'], 's': ['q'], 't': ['p'], 'u': ['t']}
LABELS = {'r': 'Root', 'v': 'V', 'o': 'O', 'n': 'N', 'm': 'M',
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
class TestInputs(unittest.TestCase):

    @test_for(check_inputs_sets)
    def test_check_inputs_sets_lst(self):
        itr, ref = check_inputs_sets(I_LST, R_LST)
        self.assertDictEqual(itr, {'a': 1, 'b': 1, 'c': 1})
        self.assertDictEqual(ref, {'a': 1, 'b': 1, 'c': 1, 'd': 1, 'e': 1,
                                   'f': 1, 'g': 1, 'h': 1, 'i': 1})

    @test_for(check_inputs_sets)
    def test_check_inputs_sets_set(self):
        itr, ref = check_inputs_sets(set(I_LST), set(R_LST))
        self.assertDictEqual(itr, {'a': 1, 'b': 1, 'c': 1})
        self.assertDictEqual(ref, {'a': 1, 'b': 1, 'c': 1, 'd': 1, 'e': 1,
                                   'f': 1, 'g': 1, 'h': 1, 'i': 1})

    @test_for(check_inputs_sets)
    def test_check_inputs_sets_dct(self):
        itr, ref = check_inputs_sets(I_DCT, R_DCT)
        self.assertDictEqual(itr, I_DCT)
        self.assertDictEqual(ref, R_DCT)

    @test_for(check_inputs_sets)
    def test_check_inputs_sets_mix(self):
        itr, ref = check_inputs_sets(I_DCT, R_LST)
        self.assertDictEqual(itr, I_DCT)
        self.assertDictEqual(ref, {'a': 1, 'b': 1, 'c': 1, 'd': 1, 'e': 1,
                                   'f': 1, 'g': 1, 'h': 1, 'i': 1})

    @test_for(check_inputs_sets)
    def test_check_inputs_sets_no_ref(self):
        itr, ref = check_inputs_sets(I_DCT, None)
        self.assertDictEqual(itr, I_DCT)
        self.assertDictEqual(ref, I_DCT)

    @test_for(check_inputs_sets)
    def test_check_inputs_sets_no_itr(self):
        self.assertRaises(ValueError, check_inputs_sets, None, None)

    @test_for(check_inputs_sets)
    def test_check_inputs_sets_wrong_types(self):
        self.assertRaises(ValueError, check_inputs_sets, 'input', None)

    @test_for(check_inputs_sets)
    def test_check_inputs_sets_wrong_types_dct(self):
        self.assertRaises(ValueError, check_inputs_sets, {'a': 8, 'b': 9, 'c': 'no'}, None)
        self.assertRaises(ValueError, check_inputs_sets, {'a': 8, 'b': 9, 'c': -2}, None)
        self.assertRaises(ValueError, check_inputs_sets, {'a': 8, 'b': 9, 'c': 0}, None)

    @test_for(get_file)
    def test_get_file_classes(self):
        tested_file = '../../ontosunburst/Inputs/ec__18jun25__classes.json'
        abs_path = os.path.abspath(tested_file)
        file = get_file('ec', 'classes.json')
        self.assertEqual(file, abs_path)

    @test_for(get_file)
    def test_get_file_labels(self):
        tested_file = '../../ontosunburst/Inputs/go_cc__22jul25__labels.json'
        abs_path = os.path.abspath(tested_file)
        file = get_file('go_cc', 'labels.json')
        self.assertEqual(file, abs_path)

    @test_for(get_file)
    def test_get_file_not_exists(self):
        self.assertRaises(FileNotFoundError, get_file, 'metacyc', 'labels.json')

    @test_for(get_ontology_dag_dict)
    def test_get_ontology_dag_dict_name(self):
        onto_dag = get_ontology_dag_dict('ec', None)
        file_dag = os.path.join('..', '..', 'ontosunburst', 'Inputs', 'ec__18jun25__classes.json')
        with open(file_dag, 'r') as f:
            expected = json.load(f)
        self.assertDictEqual(onto_dag, expected)

    @test_for(get_ontology_dag_dict)
    def test_get_ontology_dag_dict_dict(self):
        onto_dag = get_ontology_dag_dict(None, ONTO)
        self.assertDictEqual(onto_dag, ONTO)

    @test_for(get_ontology_dag_dict)
    def test_get_ontology_dag_dict_name_dict(self):
        onto_dag = get_ontology_dag_dict('ec', ONTO)
        self.assertDictEqual(onto_dag, ONTO)

    @test_for(get_ontology_dag_dict)
    def test_get_ontology_dag_dict_file(self):
        onto_file = os.path.join('test_files', 'toy_onto.json')
        onto_dag = get_ontology_dag_dict(None, onto_file)
        self.assertDictEqual(onto_dag, ONTO)

    @test_for(get_ontology_dag_dict)
    def test_get_ontology_dag_dict_wrong_name(self):
        self.assertRaises(ValueError, get_ontology_dag_dict, 'onto', None)

    @test_for(get_ontology_dag_dict)
    def test_get_ontology_dag_dict_wrong_type(self):
        self.assertRaises(ValueError, get_ontology_dag_dict, None, ['a', 'b', 'c'])

    @test_for(get_ontology_dag_dict)
    def test_get_ontology_dag_dict_not_a_file(self):
        self.assertRaises(FileNotFoundError, get_ontology_dag_dict, None, 'this_is_a_file.txt')

    @test_for(get_ontology_dag_dict)
    def test_get_ontology_dag_dict_wrong_file(self):
        onto_file = os.path.join('test_files', 'toy_onto_bad_file.txt')
        self.assertRaises(json.decoder.JSONDecodeError, get_ontology_dag_dict, None, onto_file)

    @test_for(get_ontology_dag_dict)
    def test_get_ontology_dag_dict_wrong_dict(self):
        onto_dict = {'a': ['b'], 'c': ['b', 'f'], 'd': 'a'}
        o = get_ontology_dag_dict(None, onto_dict)
        print(o)
        # self.assertRaises(json.decoder.JSONDecodeError, get_ontology_dag_dict, None, onto_file)

    @test_for(detect_cycles)
    def test_detect_cycles_ok(self):
        cycles = detect_cycles(ONTO)
        self.assertIsNone(cycles)

    @test_for(detect_cycles)
    def test_detect_cycles_cycles(self):
        onto_with_cycles = copy.deepcopy(ONTO)
        onto_with_cycles['o'].append('d')
        onto_with_cycles['x'].append('i')
        onto_with_cycles['w'].append('w')
        # Must detect ['w'], ['o', 'd', 'j'] and ['i', 'x'] cycles
        self.assertRaises(ValueError, detect_cycles, onto_with_cycles)

    @test_for(check_input_ids_to_ontology_mapping)
    def test_check_classes_concordance_ok(self):
        all_classes = set(R_LST)
        unmapped = check_input_ids_to_ontology_mapping(ONTO, all_classes)
        self.assertEqual(set(), unmapped)

    @test_for(check_input_ids_to_ontology_mapping)
    def test_check_classes_concordance_not_found(self):
        all_classes = set(R_LST)
        all_classes.add('truc')
        all_classes.add('autre truc')
        unmapped = check_input_ids_to_ontology_mapping(ONTO, all_classes)
        self.assertEqual(unmapped, {'truc', 'autre truc'})

    @test_for(get_ontology_root)
    def test_get_ontology_root_1root(self):
        root, onto_dag = get_ontology_root(ONTO)
        self.assertEqual(root, 'r')
        self.assertDictEqual(onto_dag, ONTO)

    @test_for(get_ontology_root)
    def test_get_ontology_root_several(self):
        ontology_several_roots = copy.deepcopy(ONTO)
        ontology_several_roots['n'].append('other root')
        ontology_several_roots['o'].append('another root')
        ontology_several_roots['r'] = []
        root, onto_dag = get_ontology_root(ontology_several_roots)
        self.assertEqual(root, 'Sunburst Root')
        ontology_several_roots['other root'] = ['Sunburst Root']
        ontology_several_roots['another root'] = ['Sunburst Root']
        ontology_several_roots['r'] = ['Sunburst Root']
        self.assertDictEqual(onto_dag, ontology_several_roots)
        roots = set()
        not_roots = {x for x, y in onto_dag.items() if y != []}
        for p_list in onto_dag.values():
            for p in p_list:
                if p not in not_roots:
                    roots.add(p)
        self.assertEqual(len(roots), 1)

