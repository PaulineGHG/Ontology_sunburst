import os
import json
from typing import List, Dict
from time import time
import plotly.graph_objects as go

from ontosunburst.ontology import get_abundance_dict, get_classes_abundance, get_classes_scores, \
    extract_classes, reduce_d_ontology, METACYC, CHEBI, EC, GO, KEGG, ROOTS

from ontosunburst.data_table_tree import DataTable, BINOMIAL_TEST, ROOT_CUT
from ontosunburst.sunburst_fig import generate_sunburst_fig, TOPOLOGY_A, ENRICHMENT_A

# ==================================================================================================
#                                           CONSTANTS
# ==================================================================================================

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
# Dictionary json files
# ---------------------
DEFAULT_FILE = {METACYC: os.path.join(CURRENT_DIR, 'Inputs', 'MetaCyc26_0_classes.json'),
                EC: os.path.join(CURRENT_DIR, 'Inputs', 'enzymes_ontology.json'),
                KEGG: os.path.join(CURRENT_DIR, 'Inputs', 'kegg_onto.json')}
# Names json files
# ----------------
DEFAULT_NAMES = {EC: os.path.join(CURRENT_DIR, 'Inputs', 'enzymes_class_names.json'),
                 METACYC: None, KEGG: None, CHEBI: None, GO: None}
DEFAULT = 'default'
# Sparql URL
# ----------
DEFAULT_URL = {CHEBI: 'http://localhost:3030/chebi/',
               GO: 'http://localhost:3030/go/'}


# ==================================================================================================
#                                            WORKFLOW
# ==================================================================================================

def ontosunburst(metabolic_objects: List[str],
                 ontology: str = None,
                 root: str = None,
                 abundances: List[float] = None,
                 scores: List[float] = None,
                 reference_set: List[str] = None,
                 ref_abundances: List[float] = None,
                 analysis: str = TOPOLOGY_A,
                 output: str = 'sunburst',
                 write_output: bool = True,
                 class_ontology: str or Dict[str, str] = None,
                 labels: str or Dict[str, str] = DEFAULT,
                 endpoint_url: str = None,
                 test: str = BINOMIAL_TEST,
                 root_cut: str = ROOT_CUT,
                 ref_base: bool = False,
                 show_leaves: bool = False,
                 **kwargs) -> go.Figure:
    """ Main function to be called generating the sunburst figure

    Parameters
    ----------
    metabolic_objects: List[str]
        Set of metabolic objects to classify
    ontology: str (optional, default=None)
        Ontology to use, must be in : [metacyc, ec, chebi, kegg, go]
    root: str (optional, default=None)
        Root item of the ontology.
    abundances: List[str] (optional, default=None)
        Abundance values associated to metabolic_objects list parameter
    reference_set: List[str] (optional, default=None)
        Set of reference metabolic objects
    ref_abundances: List[str] (optional, default=None)
        Abundance values associated to reference_set list parameter
    analysis: str (optional, default=topology)
        Analysis mode, must be in : [topology, enrichment]
    output: str (optional, default=None)
        Path of the output to save figure
    write_output: bool (optional, default=True)
        True to write the html figure and tsv class files, False to only return figure
    class_ontology: str or Dict[str, str] (optional, default=None)
        Class ontology dictionary or json file.
    labels: str or Dict[str, str] (optional, default=default)
        Path to ID-LABELS association json file or ID-LABELS association dictionary
    endpoint_url: str (optional, default=None)
        URL of ChEBI or GO ontology for SPARQL requests
    test: str (optional, default=binomial)
        Type of test if analysis=enrichment, must be in : [binomial, hypergeometric]
    root_cut: str (optional, default=cut)
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
    start_time = time()
    # LOAD NAMES -----------------------------------------------------------------------------------
    if labels == DEFAULT:
        if ontology is not None:
            labels = DEFAULT_NAMES[ontology]
        else:
            labels = None
    if labels is not None:
        if type(labels) == str:
            with open(labels, 'r') as f:
                names = json.load(f)
        else:
            names = labels
    else:
        names = None
    # DICTIONARY / JSON INPUT ----------------------------------------------------------------------
    if ontology == METACYC or ontology == EC or ontology == KEGG or ontology is None:
        if ontology is None:
            if class_ontology is None:
                raise ValueError('If no default ontology, must fill class_ontology parameter')
            if root is None:
                raise ValueError('If no default ontology, must fill root parameter')
        else:
            if class_ontology is None:
                class_ontology = DEFAULT_FILE[ontology]
        if type(class_ontology) == str:
            with open(class_ontology, 'r') as f:
                class_ontology = json.load(f)
    # SPARQL URL INPUT -----------------------------------------------------------------------------
    elif ontology == CHEBI or ontology == GO:
        if endpoint_url is None:
            endpoint_url = DEFAULT_URL[ontology]
    # ELSE -----------------------------------------------------------------------------------------
    else:
        raise ValueError(f'ontology parameter must be in {[METACYC, EC, KEGG, GO, CHEBI, None]}')
    # GET ROOT -------------------------------------------------------------------------------------
    if ontology is not None:
        root = ROOTS[ontology]
    # WORKFLOW -------------------------------------------------------------------------------------
    fig = _global_analysis(ontology=ontology, analysis=analysis,
                           metabolic_objects=metabolic_objects, abundances=abundances,
                           scores=scores,
                           reference_set=reference_set, ref_abundances=ref_abundances,
                           d_classes_ontology=class_ontology, endpoint_url=endpoint_url,
                           output=output, write_output=write_output, names=names,
                           test=test, root=root, root_cut=root_cut, ref_base=ref_base,
                           show_leaves=show_leaves, **kwargs)
    end_time = time()
    print(f'Execution time : {end_time - start_time} seconds')
    return fig


# ==================================================================================================
#                                             FUNCTIONS
# ==================================================================================================
def _global_analysis(ontology, analysis, metabolic_objects, abundances, scores, reference_set,
                     ref_abundances, d_classes_ontology, endpoint_url, output, write_output, names,
                     test, root, root_cut, ref_base, show_leaves, **kwargs):
    """

    Parameters
    ----------
    ontology
    analysis
    metabolic_objects
    abundances
    reference_set
    ref_abundances
    d_classes_ontology
    endpoint_url
    output
    write_output
    names
    test
    root
    root_cut
    ref_base
    show_leaves
    kwargs

    Returns
    -------

    """
    obj_all_classes, d_classes_ontology, names = extract_classes(ontology, metabolic_objects, root,
                                                                 d_classes_ontology=d_classes_ontology,
                                                                 endpoint_url=endpoint_url)
    classes_scores = None
    abundances_dict = get_abundance_dict(abundances=abundances,
                                         metabolic_objects=metabolic_objects,
                                         ref=False)
    classes_abundance = get_classes_abundance(obj_all_classes, abundances_dict, show_leaves)
    if scores is not None:
        scores_dict = get_abundance_dict(abundances=scores,
                                         metabolic_objects=metabolic_objects,
                                         ref=False)
        classes_scores = get_classes_scores(classes_abundance, scores_dict, root)

    if not obj_all_classes:
        print('No object classified, passing.')
        return None

    if write_output:
        write_met_classes(ontology, obj_all_classes, output)

    if reference_set is not None:
        ref_abundances_dict = get_abundance_dict(abundances=ref_abundances,
                                                 metabolic_objects=reference_set,
                                                 ref=True)
        ref_all_classes, d_classes_ontology, names = extract_classes(ontology, reference_set, root,
                                                              d_classes_ontology=d_classes_ontology,
                                                              endpoint_url=endpoint_url)
        ref_classes_abundance = get_classes_abundance(ref_all_classes, ref_abundances_dict,
                                                      show_leaves)
    else:
        ref_classes_abundance = None

    # Enrichment figure
    if analysis == ENRICHMENT_A:
        return _enrichment_analysis(ref_classes_abundance=ref_classes_abundance,
                                    classes_abundance=classes_abundance,
                                    classes_scores=classes_scores,
                                    d_classes_ontology=d_classes_ontology,
                                    output=output, write_output=write_output, names=names,
                                    test=test, root=root, root_cut=root_cut,
                                    ref_base=ref_base, **kwargs)

    # Proportion figure
    elif analysis == TOPOLOGY_A:
        return _topology_analysis(ref_classes_abundance=ref_classes_abundance,
                                  classes_abundance=classes_abundance,
                                  d_classes_ontology=d_classes_ontology,
                                  output=output, write_output=write_output, names=names,
                                  root=root, root_cut=root_cut, ref_base=ref_base,
                                  **kwargs)

    else:
        raise ValueError(f'Value of analysis parameter must be in : {[TOPOLOGY_A, ENRICHMENT_A]}')


def _topology_analysis(ref_classes_abundance, classes_abundance, d_classes_ontology, output,
                       write_output, names, root, root_cut, ref_base, **kwargs) -> go.Figure:
    """ Performs the topology analysis showing the abundance of each classes.

    Parameters
    ----------
    ref_classes_abundance
    classes_abundance
    d_classes_ontology
    output
    write_output
    names
    root
    root_cut
    ref_base

    Returns
    -------
    go.Figure
        Plotly graph_objects figure of the sunburst
    """
    if ref_classes_abundance is not None:
        ref_set = True
        d_classes_ontology = reduce_d_ontology(d_classes_ontology,
                                               {**ref_classes_abundance, **classes_abundance})
        data = DataTable()
        data.fill_parameters(ref_abundance=ref_classes_abundance,
                             parent_dict=d_classes_ontology,
                             root_item=root, set_abundance=classes_abundance,
                             names=names)
    else:
        ref_set = False
        d_classes_ontology = reduce_d_ontology(d_classes_ontology, classes_abundance)
        data = DataTable()
        data.fill_parameters(ref_abundance=classes_abundance,
                             parent_dict=d_classes_ontology,
                             root_item=root, names=names)
    data.calculate_proportions(ref_base)
    data.cut_root(root_cut)
    return generate_sunburst_fig(data=data, output=output, analysis=TOPOLOGY_A,
                                 root_cut=root_cut, ref_set=ref_set,
                                 write_fig=write_output, **kwargs)


def _enrichment_analysis(ref_classes_abundance, classes_abundance, classes_scores,
                         d_classes_ontology, output,
                         write_output, names, test, root, root_cut, ref_base, **kwargs) \
        -> go.Figure:
    """ Performs the enrichment analysis showing enrichment of the interest set in comparison to
    the reference set.

    Parameters
    ----------
    ref_classes_abundance
    classes_abundance
    d_classes_ontology
    output
    write_output
    names
    test
    root
    root_cut
    ref_base

    Returns
    -------
    go.Figure
        Plotly graph_objects figure of the sunburst
    """
    data = DataTable()
    if ref_classes_abundance is not None:
        d_classes_ontology = reduce_d_ontology(d_classes_ontology,
                                               {**ref_classes_abundance, **classes_abundance})
        ref_set = True
        data.fill_parameters(ref_abundance=ref_classes_abundance, parent_dict=d_classes_ontology,
                             root_item=root, set_abundance=classes_abundance, names=names,
                             ref_base=ref_base)
    else:
        ref_set = False
        d_classes_ontology = reduce_d_ontology(d_classes_ontology, classes_abundance)
        data.fill_parameters(ref_abundance=classes_abundance,
                             parent_dict=d_classes_ontology,
                             root_item=root, names=names)
    data.calculate_proportions(ref_base)
    print(data.len)
    significant = data.make_enrichment_analysis(test, classes_scores)
    print(data)
    print(significant)
    data.cut_root(root_cut)
    return generate_sunburst_fig(data=data, output=output, analysis=ENRICHMENT_A, test=test,
                                 significant=significant, ref_set=ref_set,
                                 write_fig=write_output, **kwargs)


def write_met_classes(ontology: str, all_classes: Dict[str, List[str]], output: str):
    """ Writes, for each input class, all its ancestors in a .tsv file.

    Parameters
    ----------
    ontology
    all_classes
    output
    """
    if ontology is None:
        ontology = ''
    links_dict = {METACYC: 'https://metacyc.org/compound?orgid=META&id=',
                  CHEBI: 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:',
                  EC: 'https://enzyme.expasy.org/EC/',
                  KEGG: 'https://www.genome.jp/entry/',
                  GO: 'https://amigo.geneontology.org/amigo/term/',
                  '': ''}
    with open(f'{output}.tsv', 'w') as f:
        f.write('\t'.join(['ID', 'Classes', 'Link']) + '\n')
        for met, classes, in all_classes.items():
            link = links_dict[ontology] + met
            f.write('\t'.join([met, ', '.join(classes), link]) + '\n')
