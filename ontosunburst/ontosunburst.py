import os
import json
from typing import List
import plotly.graph_objects as go

from ontosunburst.ontology import get_abundance_dict, get_classes_abondance, get_children_dict, \
    extract_classes, reduce_d_ontology, METACYC, CHEBI, EC, GO, KEGG, ROOTS

from ontosunburst.sunburst_fig import get_fig_parameters, get_data_proportion, \
    generate_sunburst_fig, BINOMIAL_TEST, TOPOLOGY_A, ENRICHMENT_A, ROOT_CUT

# ==================================================================================================
#                                           CONSTANTS
# ==================================================================================================

# DEFAULT FILES
# -------------
CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
# For MetaCyc
METACYC_FILE = os.path.join(CURRENT_DIR, 'Inputs', 'MetaCyc26_0_classes.json')
# For EC numbers
EC_ONTO_FILE = os.path.join(CURRENT_DIR, 'Inputs', 'enzymes_ontology.json')
EC_NAMES_FILE = os.path.join(CURRENT_DIR, 'Inputs', 'enzymes_class_names.json')
# For KEGG
KEGG_ONTO_FILE = os.path.join(CURRENT_DIR, 'Inputs', 'kegg_onto.json')
# For ChEBI
CHEBI_URL = 'http://localhost:3030/chebi/'
# For GO
GO_URL = 'http://localhost:3030/go/'


# ==================================================================================================
#                                            WORKFLOW
# ==================================================================================================

def ontosunburst(ontology: str,
                 metabolic_objects: List[str],
                 abundances: List[float] = None,
                 reference_set: List[str] = None,
                 ref_abundances: List[float] = None,
                 analysis: str = TOPOLOGY_A,
                 output: str = None,
                 class_file: str = None,
                 names_file: str = None,
                 endpoint_url: str = None,
                 test: str = BINOMIAL_TEST,
                 full: bool = True,
                 total: bool = True,
                 root_cut: str = ROOT_CUT,
                 ref_base: bool = False,
                 show_leaves: bool = False) -> go.Figure:
    """

    Parameters
    ----------
    ontology: str
        Ontology to use, must be in : [metacyc, ec, chebi]
    metabolic_objects: List[str]
        Set of metabolic objects to classify
    abundances: List[str] (optional, default=None)
        Abundance values associated to metabolic_objects list parameter
    reference_set: List[str] (optional, default=None)
        Set of reference metabolic objects
    ref_abundances: List[str] (optional, default=None)
        Abundance values associated to reference_set list parameter
    analysis: str (optional, default=topology)
        Analysis mode, must be in : [topology, enrichment]
    output: str (optional, default=None)
        Path to output to save figure
    class_file: str (optional, default=None)
        Path to class ontology file
    names_file: str (optional, default=None)
        Path to EC_ID - EC_NAME association json file
    endpoint_url: str (optional, default=None)
        URL of ChEBI ontology for SPARQL requests
    test: str (optional, default=binomial)
        Type of test if analysis=enrichment, must be in : [binomial, hypergeometric]
    full: bool (optional, default=True)
        True to duplicate labels if +1 parents (False to take exactly 1 random parent)
    total: bool (optional, default=True)
        True to have branch values proportional of the total parent (may not work in some cases)
    root_cut: str (optional, default=ROOT_CUT)
        mode for root cutting (uncut, cut, total)
    ref_base: bool (optional, default=False)
        True to have the base classes representation of the reference set in the figure.
    show_leaves: bool (optional, default=False)
        True to show input metabolic objets at sunburst leaves

    Returns
    -------
    go.Figure
        Plotly graph_objects figure of the sunburst
    """

    # METACYC --------------------------------------------------------------------------------------
    if ontology == METACYC:
        if class_file is None:
            class_file = METACYC_FILE
        names = None
        with open(class_file, 'r') as f:
            d_classes_ontology = json.load(f)
    # EC -------------------------------------------------------------------------------------------
    elif ontology == EC:
        if class_file is None:
            class_file = EC_ONTO_FILE
        if names_file is None:
            names_file = EC_NAMES_FILE
        with open(class_file, 'r') as f:
            d_classes_ontology = json.load(f)
        with open(names_file, 'r') as f:
            names = json.load(f)
    # KEGG -----------------------------------------------------------------------------------------
    elif ontology == KEGG:
        if class_file is None:
            class_file = KEGG_ONTO_FILE
        names = None
        with open(class_file, 'r') as f:
            d_classes_ontology = json.load(f)
    # CHEBI ----------------------------------------------------------------------------------------
    elif ontology == CHEBI:
        endpoint_url = CHEBI_URL
        d_classes_ontology = None
        names = None
    # GO -------------------------------------------------------------------------------------------
    elif ontology == GO:
        endpoint_url = GO_URL
        d_classes_ontology = None
        names = None
    # ELSE -----------------------------------------------------------------------------------------
    else:
        raise ValueError(f'ontology parameter must be in : {[METACYC, EC, KEGG, CHEBI, GO]}')
    # WORKFLOW -------------------------------------------------------------------------------------
    return _global_analysis(ontology=ontology, analysis=analysis,
                            metabolic_objects=metabolic_objects, abundances=abundances,
                            reference_set=reference_set, ref_abundances=ref_abundances,
                            d_classes_ontology=d_classes_ontology, endpoint_url=endpoint_url,
                            output=output, full=full, names=names, total=total, test=test,
                            root=ROOTS[ontology], root_cut=root_cut, ref_base=ref_base,
                            show_leaves=show_leaves)


# ==================================================================================================
#                                             FUNCTIONS
# ==================================================================================================
def _global_analysis(ontology, analysis, metabolic_objects, abundances, reference_set,
                     ref_abundances, d_classes_ontology, endpoint_url, output, full, names, total,
                     test, root, root_cut, ref_base, show_leaves):
    abundances_dict = get_abundance_dict(abundances=abundances,
                                         metabolic_objects=metabolic_objects,
                                         ref=False)

    obj_all_classes, d_classes_ontology = extract_classes(ontology, metabolic_objects, root,
                                                          d_classes_ontology=d_classes_ontology,
                                                          endpoint_url=endpoint_url)
    classes_abundance = get_classes_abondance(obj_all_classes, abundances_dict, show_leaves)

    if not obj_all_classes:
        print('No metabolic object classified, passing.')
        return None

    if output is not None:
        write_met_classes(ontology, obj_all_classes, output)

    if reference_set is not None:
        ref_abundances_dict = get_abundance_dict(abundances=ref_abundances,
                                                 metabolic_objects=reference_set,
                                                 ref=True)
        ref_all_classes, d_classes_ontology = extract_classes(ontology, reference_set, root,
                                                              d_classes_ontology=d_classes_ontology,
                                                              endpoint_url=endpoint_url)
        ref_classes_abundance = get_classes_abondance(ref_all_classes, ref_abundances_dict,
                                                      show_leaves)

    else:
        ref_classes_abundance = None

    # Enrichment figure
    if analysis == ENRICHMENT_A:
        if ref_classes_abundance is not None:
            return _enrichment_analysis(ref_classes_abundance=ref_classes_abundance,
                                        classes_abundance=classes_abundance,
                                        d_classes_ontology=d_classes_ontology,
                                        output=output, full=full, names=names, total=total,
                                        test=test, root=root, root_cut=root_cut, ref_base=ref_base,
                                        show_leaves=show_leaves)
        else:
            raise AttributeError('Missing reference set parameter')

    # Proportion figure
    elif analysis == TOPOLOGY_A:
        return _topology_analysis(ref_classes_abundance=ref_classes_abundance,
                                  classes_abundance=classes_abundance,
                                  d_classes_ontology=d_classes_ontology,
                                  output=output, full=full, names=names, total=total,
                                  root=root, root_cut=root_cut, ref_base=ref_base,
                                  show_leaves=show_leaves)

    else:
        raise ValueError(f'Value of analysis parameter must be in : {[TOPOLOGY_A, ENRICHMENT_A]}')


def _topology_analysis(ref_classes_abundance, classes_abundance, d_classes_ontology, output, full,
                       names, total, root, root_cut, ref_base, show_leaves) -> go.Figure:
    """ Performs the proportion analysis

    Parameters
    ----------
    ref_classes_abundance
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
    if ref_classes_abundance is not None:
        d_classes_ontology = reduce_d_ontology(d_classes_ontology, ref_classes_abundance)

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
        d_classes_ontology = reduce_d_ontology(d_classes_ontology, classes_abundance)
        data = get_fig_parameters(classes_abondance=classes_abundance,
                                  parent_dict=d_classes_ontology,
                                  children_dict=get_children_dict(d_classes_ontology),
                                  root_item=root, full=full, names=names)
        data = get_data_proportion(data, total)
        names = names is not None
        return generate_sunburst_fig(data=data, output=output, analysis=TOPOLOGY_A,
                                     names=names, total=total, root_cut=root_cut, ref_base=ref_base)


def _enrichment_analysis(ref_classes_abundance, classes_abundance, d_classes_ontology, output, full,
                         names, total, test, root, root_cut, ref_base, show_leaves) -> go.Figure:
    """ Performs the comparison analysis

    Parameters
    ----------
    ref_classes_abundance
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
    return generate_sunburst_fig(data=data, output=output, analysis=ENRICHMENT_A,
                                 ref_classes_abundance=ref_classes_abundance, test=test,
                                 names=names, total=total, root_cut=root_cut)


def write_met_classes(ontology, all_classes, output):
    links_dict = {METACYC: 'https://metacyc.org/compound?orgid=META&id=',
                  CHEBI: 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:',
                  EC: 'https://enzyme.expasy.org/EC/',
                  KEGG: 'https://www.genome.jp/entry/',
                  GO: 'https://amigo.geneontology.org/amigo/term/'}
    with open(f'{output}.tsv', 'w') as f:
        f.write('\t'.join(['ID', 'Classes', 'Link']) + '\n')
        for met, classes, in all_classes.items():
            link = links_dict[ontology] + met
            f.write('\t'.join([met, ', '.join(classes), link]) + '\n')
