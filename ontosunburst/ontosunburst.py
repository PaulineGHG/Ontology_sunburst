import os
import json
from typing import Collection
import plotly.graph_objects as go

from padmet.classes import PadmetRef

from ontosunburst.ontology import get_all_classes, get_classes_abondance, get_children_dict, \
    extract_chebi_roles, extract_metacyc_classes, extract_ec_classes

from ontosunburst.sunburst_fig import get_fig_parameters, get_data_proportion, \
    generate_sunburst_fig, BINOMIAL_TEST, TOPOLOGY_A, ENRICHMENT_A, ROOT_CUT

# CONSTANTS ========================================================================================

# DEFAULT FILES
# -------------
CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
# For MetaCyc
CLASS_FILE = os.path.join(CURRENT_DIR, 'Inputs/classes.json')
METACYC_FILE = os.path.join(CURRENT_DIR, 'Inputs/metacyc_26.0_prot70.padmet')
# For EC numbers
ENZYME_ONTO_FILE = os.path.join(CURRENT_DIR, 'Inputs/enzymes_ontology.json')
NAMES_FILE = os.path.join(CURRENT_DIR, 'Inputs/enzymes_class_names.json')

# ONTOLOGIES
# ----------
METACYC_ROOT = 'FRAMES'
CHEBI_ROLE_ROOT = 'role'
EC_ROOT = 'Enzyme'


# WORKFLOW =========================================================================================
def metacyc_ontosunburst(metabolic_objects: Collection[str], reference_set: Collection[str] = None,
                         analysis: str = TOPOLOGY_A,
                         output: str = None, class_file: str = CLASS_FILE,
                         padmet_ref: str = METACYC_FILE, test: str = BINOMIAL_TEST,
                         full: bool = True, total: bool = True, root_cut: str = ROOT_CUT,
                         ref_base: bool = True) -> go.Figure:
    """ Classify and plot a sunburst from a list of metabolic objects with MetaCyc ontology Ids

    Parameters
    ----------
    metabolic_objects: Collection[str]
        Set of metabolic objects to classify
    reference_set: Collection[str] (optional, default=None)
        Set of reference metabolic objects
    analysis: str (optional, default=topology)
        Analysis mode : topology or enrichment
    output: str (optional, default=None)
        Path to output to save figure
    class_file: str (optional, default=CLASS_FILE)
        Path to class ontology file
    padmet_ref: str (optional, default=METACYC_FILE)
        Path to metacyc padmet ref file
    test: str (optional, default=BINOMIAL_TEST)
        Type of test for enrichment analysis if reference_set is not None
    full: bool (optional, default=True)
        True to duplicate labels if +1 parents (False to take exactly 1 random parent)
    total: bool (optional, default=True)
        True to have branch values proportional of the total parent (may not work in some cases)
    root_cut: str (optional, default=ROOT_CUT)
        mode for root cutting (uncut, cut, total)
    ref_base: bool (optional, default=True)
        True to have the base classes representation of the reference set in the figure.

    Returns
    -------
    go.Figure
        Plotly graph_objects figure of the sunburst
    """
    # Load files
    padmet_ref = PadmetRef(padmet_ref)
    with open(class_file, 'r') as f:
        d_classes_ontology = json.load(f)

    # Extract set information
    obj_leaf_classes = extract_metacyc_classes(metabolic_objects, padmet_ref)
    obj_all_classes = get_all_classes(obj_leaf_classes, d_classes_ontology, METACYC_ROOT)
    classes_abundance = get_classes_abondance(obj_all_classes)
    if reference_set is not None:
        ref_leaf_classes = extract_metacyc_classes(reference_set, padmet_ref)
    else:
        ref_leaf_classes = None

    if output is not None:
        write_met_classes(obj_all_classes, output, padmet_ref)

    return global_analysis(analysis=analysis, ref_leaf_classes=ref_leaf_classes,
                           classes_abundance=classes_abundance,
                           d_classes_ontology=d_classes_ontology,
                           output=output, full=full, names=None, total=total, test=test,
                           root=METACYC_ROOT,
                           root_cut=root_cut, ref_base=ref_base)


def chebi_ontosunburst(chebi_ids: Collection[str], endpoint_url: str,
                       reference_set: Collection[str] = None, analysis: str = TOPOLOGY_A,
                       output: str = None, test: str = BINOMIAL_TEST, full: bool = True,
                       total: bool = True, root_cut: str = ROOT_CUT, ref_base: bool = True) \
        -> go.Figure:
    """ Classify and plot a sunburst from a list of ChEBI IDs with ChEBI roles ontology

    Parameters
    ----------
    chebi_ids: Collection[str]
        Set of ChEBI IDs to classify
    endpoint_url: str
        URL of ChEBI ontology for SPARQL requests
    reference_set: Collection[str] (optional, default=None)
        Set of reference ChEBI IDs
    analysis: str (optional, default=topology)
        Analysis mode : topology or enrichment
    output: str (optional, default=None)
        Path to output to save figure
    test: str (optional, default=BINOMIAL_TEST)
        Type of test for enrichment analysis if reference_set is not None
    full: bool (optional, default=True)
        True to duplicate labels if +1 parents (False to take exactly 1 random parent)
    total: bool (optional, default=True)
        True to have branch values proportional of the total parent (may not work in some cases)
    root_cut: str (optional, default=ROOT_CUT)
        mode for root cutting (uncut, cut, total)
    ref_base: bool (optional, default=True)
        True to have the base classes representation of the reference set in the figure.

    Returns
    -------
    go.Figure
        Plotly graph_objects figure of the sunburst
    """
    # Extract set information
    all_classes, d_roles_ontology = extract_chebi_roles(chebi_ids, endpoint_url)
    classes_abondance = get_classes_abondance(all_classes)
    if reference_set is not None:
        ref_all_classes, d_roles_ontology = extract_chebi_roles(reference_set, endpoint_url)
    else:
        ref_all_classes = None

    return global_analysis(analysis=analysis, ref_leaf_classes=ref_all_classes,
                           classes_abundance=classes_abondance,
                           d_classes_ontology=d_roles_ontology, output=output, full=full,
                           names=None, total=total, test=test, root=CHEBI_ROLE_ROOT,
                           root_cut=root_cut, ref_base=ref_base)


def ec_ontosunburst(ec_set: Collection[str], reference_set: Collection[str] = None,
                    analysis: str = TOPOLOGY_A, output: str = None,
                    class_file: str = ENZYME_ONTO_FILE, names_file: str = NAMES_FILE,
                    test: str = BINOMIAL_TEST, full: bool = True, total: bool = True,
                    root_cut: str = ROOT_CUT, ref_base: bool = True) -> go.Figure:
    """ Classify and plot a sunburst from a list of EC numbers with EC ontology Ids

    Parameters
    ----------
    ec_set: Collection[str]
        Set of EC numbers objects to classify (format "x.x.x.x" or "x.x.x.-")
    reference_set: Collection[str] (optional, default=None)
        Set of reference EC numbers
    analysis: str (optional, default=topology)
        Analysis mode : topology or enrichment
    output: str (optional, default=None)
        Path to output to save figure
    class_file: str (optional, default=ENZYME_ONTO_FILE)
        Path to class ontology file
    names_file: str (optional, default=NAMES_FILE)
        Path to EC_ID - EC_NAME association json file
    test: str (optional, default=BINOMIAL_TEST)
        Type of test for enrichment analysis if reference_set is not None
    full: bool (optional, default=True)
        True to duplicate labels if +1 parents (False to take exactly 1 random parent)
    total: bool (optional, default=True)
        True to have branch values proportional of the total parent (may not work in some cases)
    root_cut: str (optional, default=ROOT_CUT)
        mode for root cutting (uncut, cut, total)
    ref_base: bool (optional, default=True)
        True to have the base classes representation of the reference set in the figure.

    Returns
    -------
    go.Figure
        Plotly graph_objects figure of the sunburst
    """
    # Load files
    with open(class_file, 'r') as f:
        d_classes_ontology = json.load(f)
    with open(names_file, 'r') as f:
        names = json.load(f)

    # Extract set information
    ec_classes = extract_ec_classes(ec_set)
    all_classes = get_all_classes(ec_classes, d_classes_ontology, EC_ROOT)
    classes_abundance = get_classes_abondance(all_classes)
    if reference_set is not None:
        ref_leaf_classes = extract_ec_classes(reference_set)
    else:
        ref_leaf_classes = None

    return global_analysis(analysis=analysis, ref_leaf_classes=ref_leaf_classes,
                           classes_abundance=classes_abundance,
                           d_classes_ontology=d_classes_ontology, output=output, full=full,
                           names=names, total=total, test=test, root=EC_ROOT,
                           root_cut=root_cut, ref_base=ref_base)


# FUNCTIONS ========================================================================================

def global_analysis(analysis, ref_leaf_classes, classes_abundance, d_classes_ontology, output,
                    full, names, total, test, root, root_cut, ref_base):
    # Enrichment figure
    if analysis == ENRICHMENT_A:
        if ref_leaf_classes is not None:
            return enrichment_analysis(ref_leaf_classes=ref_leaf_classes,
                                       classes_abundance=classes_abundance,
                                       d_classes_ontology=d_classes_ontology,
                                       output=output, full=full, names=names, total=total,
                                       test=test, root=root, root_cut=root_cut, ref_base=ref_base)
        else:
            raise AttributeError('Missing reference set parameter')

    # Proportion figure
    elif analysis == TOPOLOGY_A:
        return topology_analysis(ref_leaf_classes=ref_leaf_classes,
                                 classes_abundance=classes_abundance,
                                 d_classes_ontology=d_classes_ontology,
                                 output=output, full=full, names=None, total=total,
                                 root=METACYC_ROOT, root_cut=root_cut, ref_base=ref_base)

    else:
        raise ValueError(f'Value of analysis parameter must be in : {[TOPOLOGY_A, ENRICHMENT_A]}')


def topology_analysis(ref_leaf_classes, classes_abundance, d_classes_ontology, output, full,
                      names, total, root, root_cut, ref_base) -> go.Figure:
    """ Performs the proportion analysis

    Parameters
    ----------
    ref_leaf_classes
    classes_abundance
    d_classes_ontology
    output
    full
    names
    total
    root
    root_cut
    ref_base

    Returns
    -------
    go.Figure
        Plotly graph_objects figure of the sunburst
    """
    if ref_leaf_classes is not None:
        if root == CHEBI_ROLE_ROOT:
            ref_all_classes = ref_leaf_classes
        else:
            ref_all_classes = get_all_classes(ref_leaf_classes, d_classes_ontology, root)
        ref_classes_abundance = get_classes_abondance(ref_all_classes)

        if ref_base:
            data = get_fig_parameters(classes_abondance=ref_classes_abundance,
                                      parent_dict=d_classes_ontology,
                                      children_dict=get_children_dict(d_classes_ontology),
                                      root_item=root, subset_abundance=classes_abundance,
                                      full=full, names=names)

        else:
            data = get_fig_parameters(classes_abondance=classes_abundance,
                                      parent_dict=d_classes_ontology,
                                      children_dict=get_children_dict(d_classes_ontology),
                                      root_item=root, full=full, names=names)

        data = get_data_proportion(data, total)
        names = names is not None
        return generate_sunburst_fig(data=data, output=output, analysis=TOPOLOGY_A,
                                     ref_classes_abundance=ref_classes_abundance, names=names,
                                     total=total, root_cut=root_cut, ref_base=ref_base)

    else:
        data = get_fig_parameters(classes_abondance=classes_abundance,
                                  parent_dict=d_classes_ontology,
                                  children_dict=get_children_dict(d_classes_ontology),
                                  root_item=root, full=full, names=names)
        data = get_data_proportion(data, total)
        names = names is not None
        return generate_sunburst_fig(data=data, output=output, analysis=TOPOLOGY_A,
                                     names=names, total=total, root_cut=root_cut, ref_base=ref_base)


def enrichment_analysis(ref_leaf_classes, classes_abundance, d_classes_ontology, output, full,
                        names, total, test, root, root_cut, ref_base: bool = True) -> go.Figure:
    """ Performs the comparison analysis

    Parameters
    ----------
    ref_leaf_classes
    classes_abundance
    d_classes_ontology
    output
    full
    names
    total
    test
    root
    root_cut
    ref_base

    Returns
    -------
    go.Figure
        Plotly graph_objects figure of the sunburst
    """
    if root == CHEBI_ROLE_ROOT:
        ref_all_classes = ref_leaf_classes
    else:
        ref_all_classes = get_all_classes(ref_leaf_classes, d_classes_ontology, root)
    ref_classes_abundance = get_classes_abondance(ref_all_classes)

    if ref_base:
        data = get_fig_parameters(classes_abondance=ref_classes_abundance,
                                  parent_dict=d_classes_ontology,
                                  children_dict=get_children_dict(d_classes_ontology),
                                  root_item=root, subset_abundance=classes_abundance,
                                  full=full, names=names)
        for k, v in data.items():
            print(k, len(v), v)
    else:
        data = get_fig_parameters(classes_abondance=classes_abundance,
                                  parent_dict=d_classes_ontology,
                                  children_dict=get_children_dict(d_classes_ontology),
                                  root_item=root, full=full, names=names)

    data = get_data_proportion(data, total)
    names = names is not None
    return generate_sunburst_fig(data=data, output=output, analysis=ENRICHMENT_A,
                                 ref_classes_abundance=ref_classes_abundance, test=test,
                                 names=names, total=total, root_cut=root_cut)


def write_met_classes(all_classes, output, pref):
    with open(f'{output}.tsv', 'w') as f:
        f.write('\t'.join(['Compound', 'Classes', 'Common names', 'MetaCyc link']) + '\n')
        for met, classes, in all_classes.items():
            try:
                name = ' / '.join(pref.dicOfNode[met].misc['COMMON-NAME'])
            except KeyError:
                name = ''
            link = f'https://metacyc.org/compound?orgid=META&id={met}'
            f.write('\t'.join([met, ', '.join(classes), name, link]) + '\n')
