from typing import List, Dict, Tuple, Set
import numpy as np
import scipy.stats as stats

# ==================================================================================================
# CONSTANTS
# ==================================================================================================

# Comparison tests
BINOMIAL_TEST = 'binomial'
HYPERGEO_TEST = 'hypergeometric'

# Analysis method
TOPOLOGY_A = 'topology'
ENRICHMENT_A = 'enrichment'

MAX_RELATIVE_NB = 1000000

# Keys
# ----
IDS = 'ID'
ONTO_ID = 'Onto ID'
PARENT = 'Parent'
LABEL = 'Label'
COUNT = 'Count'
REF_COUNT = 'Reference count'
PROP = 'Proportion'
REF_PROP = 'Reference proportion'
RELAT_PROP = 'Relative proportion'
PVAL = 'Pvalue'

# Root cut
ROOT_CUT = 'cut'
ROOT_TOTAL_CUT = 'total'
ROOT_UNCUT = 'uncut'


# ==================================================================================================
# FUNCTIONS
# ==================================================================================================

# Generate Data table
# --------------------------------------------------------------------------------------------------
def get_fig_parameters(classes_abondance: Dict[str, int], parent_dict: Dict[str, List[str]],
                       root_item, subset_abundance: Dict[str, int] = None,
                       names: Dict[str, str] = None) -> Dict[str, List]:
    """ Returns a dictionary of parameters to create the sunburst figure.

    Parameters
    ----------
    classes_abondance: Dict[str, int]
        Dictionary associating for each class the number of metabolites found belonging to the class
    parent_dict: Dict[str, List[str]]
        Dictionary associating for each class, its parents classes
    root_item: str
        Name of the root item of the ontology
    subset_abundance: Dict[str, int] = None
    names: Dict[str, str]
        Dictionary associating metabolic object ID to its Name

    Returns
    -------
    Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - onto ids : Onto ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
            - reference abundance value : Reference_count (int)
    """
    data = {IDS: list(),
            ONTO_ID: list(),
            PARENT: list(),
            LABEL: list(),
            COUNT: list(),
            REF_COUNT: list()}
    for c_onto_id, c_abundance in classes_abondance.items():
        c_sub_abundance = get_sub_abundance(subset_abundance, c_onto_id, c_abundance)
        if c_onto_id != root_item:
            if names is not None:
                try:
                    c_label = names[c_onto_id]
                except KeyError:
                    c_label = c_onto_id
            else:
                c_label = c_onto_id
            all_c_ids = get_all_ids(c_onto_id, c_onto_id, parent_dict, root_item, set())
            for c_id in all_c_ids:
                data = add_value_data(data=data,
                                      m_id=c_id,
                                      onto_id=c_onto_id,
                                      label=c_label,
                                      value=c_sub_abundance,
                                      base_value=c_abundance,
                                      parent=c_id[len(c_onto_id) + 2:])  # Remove c_label__ prefix
        else:
            data = add_value_data(data=data,
                                  m_id=c_onto_id,
                                  onto_id=c_onto_id,
                                  label=c_onto_id,
                                  value=c_sub_abundance,
                                  base_value=c_abundance,
                                  parent='')
    return data


def get_all_ids(m_id: str, n_id: str, parent_dict: Dict[str, List[str]], root: str,
                all_ids: Set[str]) -> Set[str]:
    """ Return recursively all unique IDs associated with a label. The IDs correspond to the path
    in the tree from the label to the root.

    Parameters
    ----------
    m_id: str
        Input ID
    n_id
        New ID
    parent_dict: Dict[str, List[str]]
        Dictionary associating for each class, its parents classes
    root: str
        Name of the root item of the ontology
    all_ids: Set[str]
        Set of unique IDs associated with a concept label.

    Returns
    -------
    Set[str]
        Set of all unique IDs associated with a concept label.
    """
    parents = parent_dict[m_id]
    for p in parents:
        nn_id = n_id + '__' + p
        if p == root:
            all_ids.add(nn_id)
        else:
            all_ids = get_all_ids(p, nn_id, parent_dict, root, all_ids)
    return all_ids


def get_sub_abundance(subset_abundance: Dict[str, float] or None, c_label: str,
                      c_abundance: float) -> float:
    """ Get the subset abundance of a reference concept.

    Parameters
    ----------
    subset_abundance: Dict[str, float] or None
        Dictionary associating to all the subset concepts, its abundance
    c_label: str
        Label of the reference concept
    c_abundance: float
        Abundance of the reference concept

    Returns
    -------
    Abundance of the subset concept. If reference concept not in the subset, returns numpy.nan. If
    the subset_abundance is None, returns c_abundance.
    """
    if subset_abundance is not None:
        try:
            c_sub_abundance = subset_abundance[c_label]
        except KeyError:
            c_sub_abundance = np.nan
    else:
        c_sub_abundance = c_abundance
    return c_sub_abundance


def add_value_data(data: Dict[str, List], m_id: str, onto_id: str, label: str, value: float,
                   base_value: float, parent: str) -> Dict[str, List]:
    """ Fill the data dictionary for a metabolite class.

    Parameters
    ----------
    data: Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - onto ids : Onto ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (float)
            - reference abundance value : Reference_count (float)
    m_id: str
        ID unique of the metabolite class to add
    onto_id: str
        ID in the ontology
    label: str
        Label (name) of the metabolite class to add
    value: int
        Abundance value of the metabolite class to add
    base_value: int

    parent: str
        Parent metabolite class of the metabolite class to add

    Returns
    -------
    Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - onto ids : Onto ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
            - reference abundance value : Reference_count (int)
    """
    if m_id in data[IDS]:
        raise ValueError(f'{m_id} already in data IDs, all IDs must be unique.')
    data[IDS].append(m_id)
    data[ONTO_ID].append(onto_id)
    data[LABEL].append(label)
    data[PARENT].append(parent)
    data[COUNT].append(value)
    data[REF_COUNT].append(base_value)
    return data


# Add Proportion to Data table
# --------------------------------------------------------------------------------------------------
def get_data_proportion(data: Dict[str, List], total: bool) -> Dict[str, List]:
    """ Add a proportion value for color. If total add relative proportion to +1 parent for branch
    value.

    Parameters
    ----------
    data: Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - onto ids : Onto ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
            - reference abundance value : Reference_count (int)
    total: bool
        True to have branch values proportional of the total parent

    Returns
    -------
    Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - onto ids : Onto ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
            - reference abundance value : Reference_count (int)
            - proportion : Proportion (0 < float <= 1)
            - reference proportion : Reference_proportion (0 < float <= 1)
            - branch proportion : Relative_prop
    """
    # Get total proportion
    max_abondance = int(np.nanmax(data[COUNT]))
    data[PROP] = [x / max_abondance for x in data[COUNT]]
    # Get reference proportion
    max_ref_abondance = np.max(data[REF_COUNT])
    data[REF_PROP] = [x / max_ref_abondance for x in data[REF_COUNT]]
    # Get proportion relative to +1 parent proportion for total branch value
    if total:
        data[RELAT_PROP] = [x for x in data[PROP]]
        p = ''
        data = get_relative_prop(data, p)
        # IDK WHY IT WORKS ???
        missed = [data[IDS][i] for i in range(len(data[IDS])) if data[RELAT_PROP][i] < 1]
        if missed:
            parents = {data[PARENT][data[IDS].index(m)] for m in missed}
            for p in parents:
                data = get_relative_prop(data, p)
            # missed = [data[IDS][i] for i in range(len(data[IDS])) if data[RELAT_PROP][i] < 1]
    return data


def get_relative_prop(data: Dict[str, List], p_id: str):
    """ Get recursively relative proportion of a parent children to itself. Add id to data
    Relative_proportion.

    Parameters
    ----------
    data: Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - onto ids : Onto ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
            - reference abundance value : Reference_count (int)
            - proportion : Proportion (0 < float <= 1)
            - reference proportion : Reference_proportion (0 < float <= 1)
            - branch proportion : Relative_prop
    p_id: str
        ID of the parent

    Returns
    -------
    Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - onto ids : Onto ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
            - reference abundance value : Reference_count (int)
            - proportion : Proportion (0 < float <= 1)
            - reference proportion : Reference_proportion (0 < float <= 1)
            - branch proportion : Relative_prop --> + actual children values
    """
    if p_id == '':
        prop_p = MAX_RELATIVE_NB
        count_p = max(data[REF_COUNT])
    else:
        prop_p = data[RELAT_PROP][data[IDS].index(p_id)]
        count_p = data[REF_COUNT][data[IDS].index(p_id)]
    index_p = [i for i, v in enumerate(data[PARENT]) if v == p_id]
    p_children = [data[IDS][i] for i in index_p]
    count_p_children = [data[REF_COUNT][i] for i in index_p]
    if sum(count_p_children) > count_p:
        total = sum(count_p_children)
    else:
        total = count_p
    for i, c in enumerate(p_children):
        prop_c = int((count_p_children[i] / total) * prop_p)
        data[RELAT_PROP][data[IDS].index(c)] = prop_c
    for c in p_children:
        if c in data[PARENT]:
            data = get_relative_prop(data, c)
    return data


# Enrichment analysis
# --------------------------------------------------------------------------------------------------
def get_data_enrichment_analysis(data: Dict[str, List], ref_classes_abundance: Dict[str, int],
                                 test: str) \
        -> Tuple[Dict[str, List], Dict[str, float]]:
    """ Performs statistical tests for enrichment analysis.

    Parameters
    ----------
    data: Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - onto ids : Onto ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
            - reference abundance value : Reference_count (int)
            - proportion : Proportion (0 < float <= 1)
            - reference proportion : Reference_proportion (0 < float <= 1)
            - branch proportion : Relative_prop
    ref_classes_abundance: Dict[str, int]
        Abundances of reference set classes
    test: str
        Type of test Binomial or Hypergeometric

    Returns
    -------
    Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - onto ids : Onto ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
            - reference abundance value : Reference_count (int)
            - proportion : Proportion (0 < float <= 1)
            - reference proportion : Reference_proportion (0 < float <= 1)
            - branch proportion : Relative_prop
            - p-value : P-values of enrichment analysis
    Dict[str, float]
        Dictionary of significant metabolic object label associated with their p-value
    """
    M = np.max(list(ref_classes_abundance.values()))  # M = ref set total item number

    m_list = [ref_classes_abundance[x] if x in ref_classes_abundance.keys() else 0 for x in
              data[ONTO_ID]]

    N = int(np.nanmax(data[COUNT]))  # N = interest set total item number
    n_list = data[COUNT]
    data[PVAL] = list()
    nb_classes = len(set([data[LABEL][i] for i in range(len(data[COUNT]))
                          if not np.isnan(data[COUNT][i])]))
    significant_representation = dict()
    for i in range(len(m_list)):
        if type(n_list[i]) == int:  # If count not nan (= if concept in interest set)
            # Binomial Test
            if test == BINOMIAL_TEST:
                p_val = stats.binomtest(n_list[i], N, m_list[i] / M, alternative='two-sided').pvalue
            # Hypergeometric Test
            elif test == HYPERGEO_TEST:
                p_val_upper = stats.hypergeom.sf(n_list[i] - 1, M, m_list[i], N)
                p_val_lower = stats.hypergeom.cdf(n_list[i], M, m_list[i], N)
                p_val = 2 * min(p_val_lower, p_val_upper)  # bilateral
            else:
                raise ValueError(f'test parameter must be in : {[BINOMIAL_TEST, HYPERGEO_TEST]}')
            if ((n_list[i] / N) - (m_list[i] / M)) > 0:  # If over-represented :
                data[PVAL].append(-np.log10(p_val))          # Positive log10(p-value)
            else:                                        # If under-represented :
                data[PVAL].append(np.log10(p_val))           # Negative log10(p-value)
            if p_val < 0.05 / nb_classes:  # Capture significant p-values : Bonferroni correction
                significant_representation[data[ONTO_ID][i]] = p_val.round(10)
        else:
            data[PVAL].append(np.nan)
    significant_representation = dict(
        sorted(significant_representation.items(), key=lambda item: item[1]))
    return data, significant_representation


# Topology management
# --------------------------------------------------------------------------------------------------
def data_cut_root(data: Dict[str, List], mode: str) -> Dict[str, List]:
    """ Filter data to cut (or not) the root to remove not necessary 100% represented classes.

    Parameters
    ----------
    data: Dict[str, List]
        Dictionary of figure parameters
    mode: str
        Mode of root cutting
        - uncut: doesn't cut and keep all nodes from ontology root
        - cut: keep only the lowest level 100% shared node
        - total: remove all 100% shared nodes (produces a pie at center)

    Returns
    -------
    data: Dict[str, List]
        Dictionary of figure parameters with root cut applied
    """
    if mode not in {ROOT_UNCUT, ROOT_CUT, ROOT_TOTAL_CUT}:
        raise ValueError(f'Root cutting mode {mode} unknown, '
                         f'must be in {[ROOT_UNCUT, ROOT_CUT, ROOT_TOTAL_CUT]}')
    if mode == ROOT_UNCUT:
        return data
    else:
        roots_ind = [i for i in range(len(data[IDS])) if data[RELAT_PROP][i] == MAX_RELATIVE_NB]
        roots = [data[IDS][i] for i in roots_ind]
        for root_id in roots:
            root_ind = data[IDS].index(root_id)
            for v in data.values():
                del v[root_ind]
        if mode == ROOT_TOTAL_CUT:
            data[PARENT] = ['' if x in roots else x for x in data[PARENT]]
    return data
