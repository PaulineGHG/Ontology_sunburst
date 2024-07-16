import unittest
import io
from functools import wraps

from ontosunburst.ontosunburst import *
from ontosunburst.data_table_tree import *

"""
Tests manually good file creation.
No automatic tests integrated.
"""

# ==================================================================================================
# GLOBAL
# ==================================================================================================

# Complex structure (multi-parental)
C_LST = ['a', 'b', 'c']
C_REF = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']
C_LAB = [1, 2, 3]
C_RAB = [1, 2, 3, 4, 5, 6, 7, 8]
C_ONTO = {'a': ['ab'], 'b': ['ab'], 'c': ['cde', 'cf'], 'd': ['cde'], 'e': ['cde', 'eg'],
          'f': ['cf'], 'g': ['gh', 'eg'], 'h': ['gh'],
          'ab': [ROOTS[METACYC]], 'cde': ['cdecf', 'cdeeg'], 'cf': ['cdecf'],
          'eg': [ROOTS[METACYC], 'cdeeg'], 'gh': [ROOTS[METACYC]],
          'cdecf': [ROOTS[METACYC]], 'cdeeg': ['cdeeg+'], 'cdeeg+': [ROOTS[METACYC]]}
C_LABELS = {'FRAMES': 'Root', 'cdeeg+': 'CDEEG+', 'cdeeg': 'CDEEG', 'cdecf': 'CDECF', 'gh': 'GH',
            'eg': 'EG', 'cde': 'CDE', 'cf': 'CF', 'h': 'H', 'g': 'G', 'f': 'F', 'e': 'E', 'd': 'D',
            'c': 'C', 'ab': 'AB', 'b': 'B'}

# Enrichment values
E_LST = ['02', '03', '04', '05', '08', '09']
E_LAB = [23, 20, 1, 4, 1, 1]
E_REF = ['01', '02', '03', '04', '05', '06', '07', '08', '09']
E_RAB = [14, 26, 20, 10, 20, 5, 1, 1, 3]
E_ONTO = {'01': ['00'], '02': ['00'], '03': ['00'], '04': ['00'], '05': ['01'],
          '06': ['01'], '07': ['01'], '08': ['02'], '09': ['02']}
E_LABElS = {'00': '0', '01': '1', '02': '2', '03': '3', '04': '4',
            '05': '5', '06': '6', '07': '7', '08': '8', '09': '9'}


# ==================================================================================================
# FUNCTIONS UTILS
# ==================================================================================================

def save_fig_json(fig, file):
    fig = fig.to_dict()
    with open(file, 'w') as f:
        json.dump(fig, f)


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

class TestOntosunburstCustomOnto(unittest.TestCase):

    # TOPOLOGY : CUSTOM ONTO
    @test_for(ontosunburst)
    def test_ontosunburst_1(self):
        fig = ontosunburst(metabolic_objects=C_LST, ontology=None, root='FRAMES',
                           abundances=C_LAB, reference_set=C_REF, ref_abundances=C_RAB,
                           analysis='topology', output='test1', write_output=False,
                           class_ontology=C_ONTO, labels=C_LABELS, endpoint_url=None,
                           test=BINOMIAL_TEST, total=True, root_cut=ROOT_TOTAL_CUT,
                           ref_base=True, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test1.json')
        fig = json.dumps(fig.to_dict(), sort_keys=True)
        with open(w_fig_file, 'r') as f:
            w_fig = json.dumps(json.load(f), sort_keys=True)
        self.assertEqual(fig, w_fig)

    @test_for(ontosunburst)
    def test_ontosunburst_2(self):
        fig = ontosunburst(metabolic_objects=C_REF, ontology=None, root='FRAMES',
                           abundances=C_RAB, reference_set=None, ref_abundances=None,
                           analysis='topology', output='test2', write_output=False,
                           class_ontology=C_ONTO, labels=DEFAULT, endpoint_url=None,
                           test=BINOMIAL_TEST, total=True, root_cut=ROOT_CUT,
                           ref_base=False, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test2.json')
        fig = json.dumps(fig.to_dict(), sort_keys=True)
        with open(w_fig_file, 'r') as f:
            w_fig = json.dumps(json.load(f), sort_keys=True)
        self.assertEqual(fig, w_fig)

    @test_for(ontosunburst)
    def test_ontosunburst_3(self):
        fig = ontosunburst(metabolic_objects=C_LST, ontology=None, root='FRAMES',
                           abundances=C_LAB, reference_set=C_REF, ref_abundances=C_RAB,
                           analysis='topology', output='test3', write_output=False,
                           class_ontology=C_ONTO, labels=C_LABELS, endpoint_url=None,
                           test=BINOMIAL_TEST, total=True, root_cut=ROOT_UNCUT,
                           ref_base=False, show_leaves=False, bg_color='#eeeeee')
        w_fig_file = os.path.join('test_files', 'test3.json')
        fig = json.dumps(fig.to_dict(), sort_keys=True)
        with open(w_fig_file, 'r') as f:
            w_fig = json.dumps(json.load(f), sort_keys=True)
        self.assertEqual(fig, w_fig)

    # ENRICHMENT : CUSTOM ONTO

    @test_for(ontosunburst)
    def test_ontosunburst_4(self):
        fig = ontosunburst(metabolic_objects=E_LST, ontology=None, root='00',
                           abundances=E_LAB, reference_set=E_REF, ref_abundances=E_RAB,
                           analysis='enrichment', output='test4', write_output=False,
                           class_ontology=E_ONTO, labels=E_LABElS, endpoint_url=None,
                           test=BINOMIAL_TEST, total=True, root_cut=ROOT_CUT,
                           ref_base=True, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test4.json')
        fig = json.dumps(fig.to_dict(), sort_keys=True)
        with open(w_fig_file, 'r') as f:
            w_fig = json.dumps(json.load(f), sort_keys=True)
        self.assertEqual(fig, w_fig)

    @test_for(ontosunburst)
    def test_ontosunburst_5(self):
        fig = ontosunburst(metabolic_objects=E_LST, ontology=None, root='00',
                           abundances=E_LAB, reference_set=E_REF, ref_abundances=E_RAB,
                           analysis='enrichment', output='test5', write_output=True,
                           class_ontology=E_ONTO, labels=E_LABElS, endpoint_url=None,
                           test=HYPERGEO_TEST, total=False, root_cut=ROOT_UNCUT,
                           ref_base=False, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test5.json')
        fig = json.dumps(fig.to_dict(), sort_keys=True)
        with open(w_fig_file, 'r') as f:
            w_fig = json.dumps(json.load(f), sort_keys=True)
        self.assertEqual(fig, w_fig)

    @test_for(ontosunburst)
    def test_ontosunburst_5(self):
        fig = ontosunburst(metabolic_objects=E_LST, ontology=None, root='00',
                           abundances=E_LAB, reference_set=E_REF, ref_abundances=E_RAB,
                           analysis='enrichment', output='test5', write_output=False,
                           class_ontology=E_ONTO, labels=E_LABElS, endpoint_url=None,
                           test=HYPERGEO_TEST, total=False, root_cut=ROOT_UNCUT,
                           ref_base=False, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test5.json')
        fig = json.dumps(fig.to_dict(), sort_keys=True)
        with open(w_fig_file, 'r') as f:
            w_fig = json.dumps(json.load(f), sort_keys=True)
        self.assertEqual(fig, w_fig)


MET_LST = ['CPD-24674', 'CPD-24687', 'CPD-24688']
MET_REF = ['CPD-24674', 'CPD-24687', 'CPD-24688',
           'CPD-12782', 'CPD-12784', 'CPD-12787',
           'CPD-12788', 'CPD-12789', 'CPD-12796',
           'CPD-12797', 'CPD-12798', 'CPD-12805',
           'CPD-12806', 'CPD-12812', 'CPD-12816',
           'CPD-1282', 'CPD-12824', 'CPD-1283']
MET_LAB = [1, 2, 3]
MET_RAB = [1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6, 7, 8, 9]

RXN_LST = ['CROTCOALIG-RXN', 'CYSTHIOCYS-RXN', 'NQOR-RXN']
RXN_REF = ['CROTCOALIG-RXN', 'CYSTHIOCYS-RXN', 'NQOR-RXN',
           'RXN-14859', 'RXN-14873', 'RXN-14920',
           'RXN-14939', 'RXN-14975', 'RXN-21632',
           'RXN-21638', 'RXN-21652', 'RXN-8954']

PWY_LST = ['2ASDEG-PWY', '4AMINOBUTMETAB-PWY', 'ALLANTOINDEG-PWY']
PWY_REF = ['2ASDEG-PWY', '4AMINOBUTMETAB-PWY', 'ALLANTOINDEG-PWY',
           'CRNFORCAT-PWY', 'PWY-7195', 'PWY-7219',
           'PWY-7251', 'PWY-7351', 'PWY-7401',
           'PWY18C3-22', 'PWY0-1600', 'SERDEG-PWY']


class TestOntosunburstMetaCyc(unittest.TestCase):

    @test_for(ontosunburst)
    def test_ontosunburst_mc1(self):
        fig = ontosunburst(metabolic_objects=MET_LST, ontology=METACYC,
                           abundances=MET_LAB, reference_set=MET_REF, ref_abundances=MET_RAB,
                           analysis=TOPOLOGY_A, output='mc1', write_output=False,
                           class_ontology=None, labels=DEFAULT, endpoint_url=None,
                           test=HYPERGEO_TEST, total=True, root_cut=ROOT_CUT,
                           ref_base=True, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_mc1.json')
        fig = json.dumps(fig.to_dict(), sort_keys=True)
        with open(w_fig_file, 'r') as f:
            w_fig = json.dumps(json.load(f), sort_keys=True)
        self.assertEqual(fig, w_fig)

    @test_for(ontosunburst)
    def test_ontosunburst_mc2(self):
        fig = ontosunburst(metabolic_objects=MET_LST, ontology=METACYC,
                           abundances=MET_LAB, reference_set=MET_REF, ref_abundances=MET_RAB,
                           analysis=TOPOLOGY_A, output='mc2', write_output=False,
                           class_ontology=None, labels={'Acids': 'ACIDS !!'}, endpoint_url=None,
                           test=HYPERGEO_TEST, total=True, root_cut=ROOT_CUT,
                           ref_base=False, show_leaves=False)
        w_fig_file = os.path.join('test_files', 'test_mc2.json')
        fig = json.dumps(fig.to_dict(), sort_keys=True)
        with open(w_fig_file, 'r') as f:
            w_fig = json.dumps(json.load(f), sort_keys=True)
        self.assertEqual(fig, w_fig)

    @test_for(ontosunburst)
    def test_ontosunburst_mc3(self):
        fig = ontosunburst(metabolic_objects=MET_REF, ontology=METACYC,
                           abundances=None, reference_set=None, ref_abundances=None,
                           analysis=TOPOLOGY_A, output='mc3', write_output=False,
                           class_ontology=None, labels=DEFAULT, endpoint_url=None,
                           test=HYPERGEO_TEST, total=True, root_cut=ROOT_CUT,
                           ref_base=False, show_leaves=True)
        w_fig_file = os.path.join('test_files', 'test_mc3.json')
        fig = json.dumps(fig.to_dict(), sort_keys=True)
        with open(w_fig_file, 'r') as f:
            w_fig = json.dumps(json.load(f), sort_keys=True)
        self.assertEqual(fig, w_fig)

    @test_for(ontosunburst)
    def test_ontosunburst_mc4(self):
        fig = ontosunburst(metabolic_objects=MET_LST, ontology=METACYC,
                           abundances=MET_LAB, reference_set=MET_REF, ref_abundances=MET_RAB,
                           analysis=ENRICHMENT_A, output='mc4', write_output=True,
                           class_ontology=None, labels=DEFAULT, endpoint_url=None,
                           test=HYPERGEO_TEST, total=True, root_cut=ROOT_CUT,
                           ref_base=False, show_leaves=False)
        w_fig_file = os.path.join('test_files', 'test_mc4.json')
        save_fig_json(fig, w_fig_file)
        fig = json.dumps(fig.to_dict(), sort_keys=True)
        with open(w_fig_file, 'r') as f:
            w_fig = json.dumps(json.load(f), sort_keys=True)
        self.assertEqual(fig, w_fig)