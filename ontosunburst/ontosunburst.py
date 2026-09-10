import copy
import json
import logging
import os
import networkx
from time import time
from typing import Tuple, TypeAlias, Literal, cast, get_args

import plotly.graph_objects as go

from ontosunburst.dag2tree import *
from ontosunburst.onto2dag import *
from ontosunburst.tree2sunburst import generate_sunburst_fig, TOPOLOGY_A, ENRICHMENT_A

# ==================================================================================================
#                                             TYPES
# ==================================================================================================
Input: TypeAlias = List[str] | Set[str] | Dict[str, float]
OntologyName: TypeAlias = Literal['metacyc', 'ec', 'chebi', 'chebi_r', 'go_cc', 'go_mf', 'go_bp',
                                  'go', 'kegg']
FileSuffix: TypeAlias = Literal['classes.json', 'labels.json']
OntologyDAG: TypeAlias = Dict[str, List[str]]
FilePath: TypeAlias = str


# ==================================================================================================
#                                           CONSTANTS
# ==================================================================================================
CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
DEFAULT_PATH = os.path.join(CURRENT_DIR, 'Inputs')

METACYC, EC, CHEBI, CHEBI_R, GO_CC, GO_MF, GO_BP, GO, KEGG = get_args(OntologyName)
ROOTS = {METACYC: 'FRAMES',
         CHEBI: 'chebi',
         CHEBI_R: 'CHEBI:50906',
         EC: 'Enzyme',
         GO_CC: 'GO:0005575',
         GO_BP: 'GO:0008150',
         GO_MF: 'GO:0003674',
         GO: 'GO',
         KEGG: 'kegg'}
CLASSES_SUFFIX, LABELS_SUFFIX = get_args(FileSuffix)

logging.basicConfig(level=logging.INFO)


# ==================================================================================================
#                                            WORKFLOW
# ==================================================================================================
def ontosunburst(interest: Input,
                 reference: Input | None = None,
                 ontology: OntologyName = None,
                 analysis: str = TOPOLOGY_A,
                 output: FilePath = 'sunburst',
                 scores: Dict[str, float] = None,
                 write_output: bool = True,
                 ontology_dag_input: str | OntologyDAG = None,
                 id_to_label_input: str or Dict[str, str] = None,
                 use_labels: bool = True,
                 test: str = BINOMIAL_TEST,
                 root_cut: str = ROOT_CUT,
                 path_cut: str = PATH_UNCUT,
                 ref_base: bool = False,
                 hide_leaves: bool = False,
                 **kwargs) -> go.Figure:
    """ Main function to be called generating the sunburst figure

    Parameters
    ----------
    interest:  list[str] | set[str] | dict[str, float]
        Interest list or set of concepts IDs to classify.
        Can be associated to a weigh as a dictionary.
    reference: list[str] | set[str] | dict[str, float] | None (optional, default=None)
        Reference list or set of concepts IDs to classify.
        Can be associated to a weigh as a dictionary.
    ontology: str (optional, default=None, values in ['metacyc', 'ec', 'chebi', 'chebi_r', 'kegg',
                                                      'go_cc', 'go_bp', 'go_mf', 'go', None])
        Ontology name to use.
    analysis: str (optional, default='topology', values in ['topology', 'enrichment'])
        Analysis mode : topology or enrichment.
    output: str (optional, default='sunburst')
        Path of the output to save figure, if None, outputs will be sunburst.html and sunburst.tsv
        files
    scores: Dict[str, float] (optional, default=None)
        Dictionary associating for each ontology ID, its precalculated enrichment score. If None
        enrichment will be calculated.
    write_output: bool (optional, default=True)
        True to write the html figure and tsv class files, False to only return plotly sunburst
        figure.
    ontology_dag_input: str or Dict[str, str] (optional, default=None)
        Ontology DAG dictionary or json file. Use if tailored ontology or alternative
        (modified, updated, ...) default ontology DAG.
    id_to_label_input: str or Dict[str, str] (optional, default=None)
        Path to ID-LABELS association json file or ID-LABELS association dictionary.
        If None default files will be used. Use if tailored ontology or alternative
        (modified, updated, ...) default ontology.
    labels: bool (optional, default=True)
        True to show labels as sunburst sectors labels, False to show ID as sunburst sectors labels.
    test: str (optional, default='binomial', values in ['binomial', 'hypergeometric'])
        Type of test if analysis=enrichment, binomial or hypergeometric test.
    root_cut: str (optional, default='cut', values in ['uncut', 'cut', 'total'])
        mode for root cutting (uncut, cut or total)
    path_cut: str (optional, default='uncut', values in ['uncut', 'deeper', 'higher', 'bound'])
        mode for nested path cutting (uncut, deeper, higher or bound)
    ref_base: bool (optional, default=False)
        True to have the base classes representation of the reference set in the figure.
    show_leaves: bool (optional, default=False)
        True to show input metabolic objets at sunburst leaves
    **kwargs

    Returns
    -------
    go.Figure
        Plotly graph_objects figure of the sunburst
    """
    start_time = time()
    # MANAGE INPUTS
    interest, reference = check_inputs_sets(interest, reference)
    all_concepts = set(interest.keys()).union(set(reference.keys()))
    # LOAD ONTOLOGY DAG DICTIONARY -----------------------------------------------------------------
    ontology_dag = get_ontology_dag_dict(ontology, ontology_dag_input)
    detect_cycles(ontology_dag)
    check_input_ids_to_ontology_mapping(ontology_dag, all_concepts)
    # GET ROOT -------------------------------------------------------------------------------------
    root, ontology_dag = get_ontology_root(ontology_dag)
    # LOAD ID TO LABELS DICTIONARY -----------------------------------------------------------------
    id_to_label = get_id_to_label_dict(id_to_label_input, ontology)


    # WORKFLOW -------------------------------------------------------------------------------------

    end_time = time()
    print(f'Execution time : {end_time - start_time} seconds')
    # return fig


def _global_analysis(analysis, interest_concepts, abundances, scores, reference_concepts,
                     ref_abundances, ontology_dag, output, write_output, id_to_label,
                     test, root, root_cut, path_cut, ref_base, show_leaves, **kwargs):
    """

    Parameters
    ----------
    analysis
    interest_concepts
    abundances
    reference_concepts
    ref_abundances
    ontology_dag
    output
    write_output
    id_to_label
    test
    root
    root_cut
    path_cut
    ref_base
    show_leaves
    kwargs

    Returns
    -------

    """
    # ONTOLOGY TO WEIGHTED DAG
    # =============================================================================================
    # Calculate all concepts weights --------------------------------------------------------------
    dag = StdDAG(concepts=interest_concepts, abundances=abundances,
                 ref_concepts=reference_concepts, ref_abundances=ref_abundances,
                 ontology_dag=ontology_dag, root=root, labels=id_to_label)

    calculated_weights = ontology_to_weighted_dag(concepts=interest_concepts, abundances=abundances,
                                                  root=root, ontology_dag=ontology_dag,
                                                  show_lvs=show_leaves)

    ref_set = reference_concepts is not None
    if ref_set:
        ref_calculated_weights = ontology_to_weighted_dag(concepts=reference_concepts,
                                                          abundances=ref_abundances, root=root,
                                                          ontology_dag=ontology_dag,
                                                          show_lvs=show_leaves)
    else:
        ref_calculated_weights = calculated_weights

    # Scores
    if ref_base:
        classes_scores = get_classes_scores(ref_calculated_weights, scores, root)
    else:
        classes_scores = get_classes_scores(calculated_weights, scores, root)

    # Reduce ontology (get DAG subgraph) ----------------------------------------------------------
    if ref_base:
        ontology_dag = reduce_d_ontology(ontology_dag, ref_calculated_weights)
        id_to_label = reduce_d_ontology(id_to_label, ref_calculated_weights)
    else:
        ontology_dag = reduce_d_ontology(ontology_dag, calculated_weights)
        id_to_label = reduce_d_ontology(id_to_label, calculated_weights)

    # if write_output:
    #     write_concepts_classes(ontology, concepts_all_classes, output, id_to_label)

    # DAG TO TREE
    # =============================================================================================
    tree_data = TreeData()
    tree_data.dag_to_tree(set_abundance=calculated_weights, ref_abundance=ref_calculated_weights,
                          parent_dict=ontology_dag, root_item=root, names=id_to_label,
                          ref_base=ref_base)

    tree_data.calculate_proportions(ref_base)
    significant = None
    if analysis == ENRICHMENT_A:
        significant = tree_data.make_enrichment_analysis(test, classes_scores)
    tree_data.cut_root(root_cut)
    tree_data.cut_nested_path(path_cut, ref_base)

    # TREE TO SUNBURST
    # =============================================================================================
    return generate_sunburst_fig(data=tree_data, output=output, analysis=analysis, test=test,
                                 significant=significant, ref_set=ref_set,
                                 write_fig=write_output, **kwargs)


# ==================================================================================================
#                                             FUNCTIONS
# ==================================================================================================
# Interest and reference inputs management
def check_valid_literal(item: str, lit):
    if item not in get_args(lit):
        raise ValueError(f'Invalid {item} argument. Must be in : {get_args(lit)}')


def check_inputs_sets(interest: Input, reference: Input | None) -> Tuple[Input, Input]:
    """ Checks for valid interest and reference parameters and raises ValueError if not.
    Converts to dictionary associating each ID to its weight (set to 1 by default).
    If no reference given, reference is equal to  interest.

    Parameters
    ----------
    interest: list[str] | set[str] | dict[str, float]
    reference: list[str] | set[str] | dict[str, float] | None

    Returns
    -------
    tuple[list[str] | set[str] | dict[str, float], list[str] | set[str] | dict[str, float]]
    """
    if interest is None:
        logging.critical('No interest set given in "interest" field.')
        raise ValueError('No interest set given in "interest" field.')
    elif type(interest) == list or type(interest) == set:
        logging.info('No abundances for interest set, '
                     'default value "1" will be used for each concept.')
        interest = {str(x): 1 for x in interest}
    elif type(interest) == dict:
        # Test weights values
        for c, w in interest.items():
            if (type(w) != float and type(w) != int) or w <= 0:
                logging.error(f'Invalid weight "{w}" for "{c}" concept in interest. '
                              f'All weights must be strictly positive numerals.')
                raise ValueError(f'Invalid weight "{w}" for "{c}" concept in interest. '
                                 f'All weights must be strictly positive numerals.')
        interest = {str(x): y for x, y in interest.items()}
    else:
        logging.error('No valid interest set given.')
        raise ValueError('No valid interest set given.')

    # Reference
    if reference is None:
        logging.info('Running ontosunburst with no reference.')
        reference = copy.deepcopy(interest)  # Maybe no need deep copy
    elif type(reference) == list or type(reference) == set:
        logging.info('No abundances for reference set, '
                     'default value "1" will be used for each concept.')
        reference = {str(x): 1 for x in reference}
    elif type(reference) == dict:
        # Test weights values
        for c, w in reference.items():
            if (type(w) != float and type(w) != int) or w <= 0:
                logging.error(f'Invalid weight "{w}" for "{c}" concept in reference. '
                              f'All weights must be strictly positive numerals.')
                raise ValueError(f'Invalid weight "{w}" for "{c}" concept in reference. '
                                 f'All weights must be strictly positive numerals.')
        reference = {str(x): y for x, y in reference.items()}
    else:
        logging.error('No valid reference set given.')
        raise ValueError('No valid reference set given.')

    return interest, reference


# Ontologies management
def get_file(ontology: OntologyName, suffix: FileSuffix) -> FilePath:
    """ Return the implemented ontology or labels file path from the ontology name and the suffix
    (classes or labels) with the current version implemented.

    Parameters
    ----------
    ontology: str
        Name of the ontology (chebi, ec, metacyc, kegg, go, go_cc, go_mf, go_bp)
    suffix: str
        Suffix of the file 'classes.json' or 'labels.json'

    Returns
    -------
    str
        File path
    """
    for file in os.listdir(DEFAULT_PATH):
        if file.startswith(ontology + '__') and file.endswith('__' + suffix):
            return os.path.join(DEFAULT_PATH, file)
    expected_file = f'{ontology}__[version]__{suffix}'
    logging.error(f'Cannot find {expected_file} file like in {DEFAULT_PATH} path.')
    raise FileNotFoundError(f'Cannot find {expected_file} file like in {DEFAULT_PATH} path.')


def merge_go_ontologies(suffix: FileSuffix) -> OntologyDAG:
    """ Merge all the 3 GO ontologies (go_cc, go_mf, go_bp) to one ontology "go" linked by a new
    root "GO".

    Parameters
    ----------
    suffix: str
        Suffix of go files to aggregate either the GO ontology classes files or labels file

    Returns
    -------
    dict[str, list[str]]
        GO Ontology DAG containing all the 3 GO ontologies merged.
    """
    go_aggregated = dict()
    for sub_go_ontology in [GO_BP, GO_CC, GO_MF]:
        dict_sub_onto_input = get_file(cast(OntologyName, sub_go_ontology), suffix)
        with open(dict_sub_onto_input, 'r') as f:
            dict_sub_onto = json.load(f)
        go_aggregated.update(dict_sub_onto)
    if suffix == CLASSES_SUFFIX:
        for sub_go_ontology in [GO_BP, GO_CC, GO_MF]:
            go_aggregated[ROOTS[sub_go_ontology]] = [ROOTS[GO]]
    return go_aggregated


def get_ontology_dag_dict(ontology: OntologyName | None,
                          ontology_dag_input: OntologyDAG | str | None) -> OntologyDAG:
    """ Return the ontology DAG as a dict depending on the inputs given (ontology name, custom
    ontology dict or custom ontology json file path)

    Parameters
    ----------
    ontology: str | None
        A default ontology name (chebi, metacyc, ec, ...) or None for custom ontology
    ontology_dag_input: dict[str, list[str]] | str | None
        A custom ontology as a dict or a json file path or None for default ontology

    Returns
    -------
    dict[str, list[str]]
        The ontology DAG as a dict.
    """
    # Case ontology_dag_input parameter not filled (default : None)
    if ontology_dag_input is None:
        # Case no default ontology : raises an error
        if ontology is None:
            logging.error('If no default ontology, must fill ontology_dag_input parameter')
            raise ValueError('If no default ontology, must fill ontology_dag_input parameter')
        # Case default ontology : get default ontology file path
        else:
            check_valid_literal(ontology, OntologyName)
            if ontology == GO:
                return merge_go_ontologies(CLASSES_SUFFIX)
            ontology_dag_input = get_file(ontology, CLASSES_SUFFIX)
            logging.info(f'Using {ontology_dag_input} file as ontology DAG.')
    # Case ontology_dag_input parameter is a file path (str)
    if type(ontology_dag_input) == str:
        if os.path.exists(ontology_dag_input):
            with open(ontology_dag_input, 'r') as f:
                ontology_dag = json.load(f)
                return ontology_dag
        else:
            logging.error(f'No file {ontology_dag_input} found, '
                          f'check for valid ontology_dag_input parameter given.')
            raise FileNotFoundError(f'No file {ontology_dag_input} found, '
                                    f'check for valid ontology_dag_input parameter given.')
    # Case ontology_dag_input parameter is a dictionary (dict)
    elif type(ontology_dag_input) == dict:
        return ontology_dag_input
    # Case ontology_dag_input parameter is not a dictionary (dict), neither a file path (str) :
    # raises an error
    else:
        raise ValueError('ontology_dag_input parameter must be a json file path (str) or a '
                         'dictionary')


def detect_cycles(onto_dag: OntologyDAG):
    """ Checks if an ontology DAG is actually a DAG by detecting cycles. Il cycles are detected,
    logs the cycles to remove and raises an error.

    Parameters
    ----------
    onto_dag: dict[str, list[str]]
        Ontology DAG dictionary
    """
    graph = networkx.DiGraph(onto_dag)
    if not networkx.is_directed_acyclic_graph(graph):
        cycles = networkx.simple_cycles(graph)
        logging.error('Cycles detected in ontology graph, cannot be used.')
        for cycle in cycles:
            logging.error(f'Cycle: {cycle} detected : fix it by removing a relation.')
        raise ValueError(f'Ontology graph given is not a DAG. Remove cycles to use.')
    else:
        logging.info('No cycles detected in ontology graph.')


def check_input_ids_to_ontology_mapping(onto_dag: OntologyDAG, all_classes: Set[str]) -> Set[str]:
    """ Checks the coverage of inputs IDs mappable on the Ontology DAG used. Warns about unmapped
    concepts. Returns the set of the unmapped concepts.

    Parameters
    ----------
    onto_dag: dict[str, list[str]]
        Ontology DAG
    all_classes: set[str]
        Set of all concepts ids (from interest and reference)

    Returns
    -------
    set[str]
        Set of IDs unmapped on the ontology DAG.
    """
    nb_concepts = len(all_classes)
    logging.info(f'{nb_concepts} concepts to map.')
    mapped = 0
    unmapped = set()
    for c in all_classes:
        if c in onto_dag:
            mapped += 1
        else:
            logging.warning(f'Concept "{c}" not found in ontology DAG.')
            unmapped.add(c)
    logging.info(f'{mapped}/{nb_concepts} concepts mapped in ontology.')
    return unmapped


def get_ontology_root(onto_dag: OntologyDAG) -> Tuple[str, Dict[str, List[str]]]:
    """ Gets the unique root concept of an ontology DAG (concept with no parent concept). If
    several roots found, creates a new root over the others. Returns the unique root and the
    modified ontology DAG if a new root has been created.

    Parameters
    ----------
    onto_dag: dict[str, list[str]]
        Ontology DAG

    Returns
    -------
    tuple[str, dict[str, list[str]]]
        Unique root ID and the Ontology DAG (eventually modified)
    """
    roots = set()
    not_roots = {x for x, y in onto_dag.items() if y != []}
    for p_list in onto_dag.values():
        for p in p_list:
            if p not in not_roots:
                roots.add(p)
    if len(roots) == 1:
        unique_root = roots.pop()
        logging.info(f'Using "{unique_root}" as ontology DAG root.')
    else:
        unique_root = 'Sunburst Root'
        logging.warning(f'Several roots found in ontology DAG. Roots: {roots}')
        logging.warning(f'Creating new unique root "{unique_root}" over the roots.')
        for root in roots:
            onto_dag[root] = [unique_root]
    return unique_root, onto_dag


def get_id_to_label_dict(id_to_label_input: Dict[str, str] | str | None,
                         ontology: OntologyName | None):
    # Case default ontology AND use of default labels file
    if ontology is not None and id_to_label_input is None:
        if ontology == GO:
            return merge_go_ontologies(LABELS_SUFFIX)
        id_to_label_input = get_file(ontology, LABELS_SUFFIX)
    # Case id_to_label_input parameter filled
    if id_to_label_input is not None:
        # Case id_to_label_input parameter is a file path (str)
        if type(id_to_label_input) == str:
            with open(id_to_label_input, 'r') as f:
                id_to_label = json.load(f)
                return id_to_label
        # Case id_to_label_input parameter is a dictionary (dict)
        elif type(id_to_label_input) == dict:
            return id_to_label_input
        # Case id_to_label_input parameter is not a dictionary (dict), neither a file
        # path (str) : raises an error
        else:
            raise ValueError('id_to_label_input parameter must be a json file path (str) or a '
                             'dictionary')


def write_concepts_classes(ontology: str, all_classes: Dict[str, Set[str]], output: str,
                           id_to_label: Dict[str, str]):
    """ Writes, for each input class, all its ancestors in a .tsv file.

    Parameters
    ----------
    ontology
    all_classes
    output
    id_to_label
    """
    if ontology is None:
        ontology = ''
    links_dict = {METACYC: 'https://metacyc.org/compound?orgid=META&id=',
                  CHEBI: 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:',
                  CHEBI_R: 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:',
                  EC: 'https://enzyme.expasy.org/EC/',
                  KEGG: 'https://www.genome.jp/entry/',
                  GO_MF: 'https://amigo.geneontology.org/amigo/term/',
                  GO_CC: 'https://amigo.geneontology.org/amigo/term/',
                  GO_BP: 'https://amigo.geneontology.org/amigo/term/',
                  GO: 'https://amigo.geneontology.org/amigo/term/',
                  '': ''}
    with open(f'{output}.tsv', 'w') as f:
        f.write('\t'.join(['ID', 'Label', 'Classes ID', 'Classes Label', 'Link']) + '\n')
        for met_id, classes_id, in all_classes.items():
            link = links_dict[ontology] + met_id
            met_lab = get_name(met_id, id_to_label)
            classes_lab = [get_name(cl, id_to_label) for cl in classes_id]
            f.write('\t'.join([met_id, met_lab, ', '.join(classes_id), ', '.join(classes_lab),
                               link]) + '\n')