import unittest
from unittest.mock import patch
import io
import sys
from ontosunburst.ontology import *

"""
Tests manually good file creation.
No automatic tests integrated.
"""


class DualWriter(io.StringIO):
    def __init__(self, original_stdout):
        super().__init__()
        self.original_stdout = original_stdout

    def write(self, s):
        super().write(s)
        self.original_stdout.write(s)


# METACYC
# ==================================================================================================

MET_LST = ['a', 'b', 'c']
MET_REF = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
MET_LAB = [1, 2, 3]
MET_RAB = [1, 2, 3, 4, 5, 6, 7, 8]
MC_ONTO = {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf'], 'd': ['cde'], 'e': ['cde', 'eg'],
           'f': ['cf'], 'g': ['gh', 'eg'], 'h': ['gh'],
           'ab': [ROOTS[METACYC]], 'cde': ['cdecf', 'cdeeg'], 'cf': ['cdecf'],
           'eg': [ROOTS[METACYC], 'cdeeg'], 'gh': [ROOTS[METACYC]],
           'cdecf': [ROOTS[METACYC]], 'cdeeg': ['cdeeg+'], 'cdeeg+': [ROOTS[METACYC]]}


# EC
# ==================================================================================================
EC_LST = ['1.4.5.6', '1.4.6.7', '2.1.2.3', '1.5.3', '1.6.9.-', '1.-.-.-', '1.4.-.-']
EC_ONTO = {'1.4.5.-': ['1.4.-.-'], '1.4.6.-': ['1.4.-.-'], '2.1.2.-': ['2.1.-.-'],
           '1.5.3.-': ['1.5.-.-'], '1.6.9.-': ['1.6.-.-'],
           '1.4.-.-': ['1.-.-.-'], '2.1.-.-': ['2.-.-.-'], '1.5.-.-': ['1.-.-.-'],
           '1.6.-.-': ['1.-.-.-'], '1.-.-.-': [ROOTS[EC]], '2.-.-.-': [ROOTS[EC]]}


class Test(unittest.TestCase):
    # CLASS EXTRACTION
    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_metacyc_classes_1(self, mock_stdout):
        d_obj = extract_metacyc_classes(MET_LST, MC_ONTO)
        output = mock_stdout.getvalue().strip()
        self.assertEqual(output, '3 metabolic objects to classify\n'
                                 '3/3 metabolic objects classified')
        self.assertEqual(d_obj, {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf']})

    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_metacyc_classes_2(self, mock_stdout):
        d_obj = extract_metacyc_classes(MET_LST + ['x'], MC_ONTO)
        output = mock_stdout.getvalue().strip()
        self.assertEqual(output, '4 metabolic objects to classify\n'
                                 'x not classified.\n'
                                 '3/4 metabolic objects classified')
        self.assertEqual(d_obj, {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf']})

    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_ec_classes_1(self, mock_stdout):
        d_obj, d_onto = extract_ec_classes(EC_LST, EC_ONTO)
        output = mock_stdout.getvalue().strip()
        wanted_d_obj = {'1.4.5.6': ['1.4.5.-'], '1.4.6.7': ['1.4.6.-'], '2.1.2.3': ['2.1.2.-'],
                        '1.5.3': ['1.5.-.-'], '1.6.9.-': ['1.6.-.-'], '1.-.-.-': ['Enzyme'],
                        '1.4.-.-': ['1.-.-.-']}
        self.assertEqual(output, '7 EC numbers to classify\n'
                                 '7/7 EC numbers classified')
        self.assertEqual(d_obj, wanted_d_obj)
        self.assertEqual(d_onto, {**EC_ONTO, **wanted_d_obj})

    @patch('sys.stdout', new_callable=lambda: DualWriter(sys.stdout))
    def test_extract_ec_classes_2(self, mock_stdout):
        d_obj, d_onto = extract_ec_classes(EC_LST + ['3.5.6.9', 'ecID'], EC_ONTO)
        output = mock_stdout.getvalue().strip()
        wanted_d_obj = {'1.4.5.6': ['1.4.5.-'], '1.4.6.7': ['1.4.6.-'], '2.1.2.3': ['2.1.2.-'],
                        '1.5.3': ['1.5.-.-'], '1.6.9.-': ['1.6.-.-'], '1.-.-.-': ['Enzyme'],
                        '1.4.-.-': ['1.-.-.-']}
        print(d_obj)
        self.assertEqual(output, '9 EC numbers to classify\n'
                                 '3.5.6.9 not classified\n'
                                 'ecID not classified\n'
                                 '7/9 EC numbers classified')
        self.assertEqual(d_obj, wanted_d_obj)
        self.assertEqual(d_onto, {**EC_ONTO, **wanted_d_obj})
