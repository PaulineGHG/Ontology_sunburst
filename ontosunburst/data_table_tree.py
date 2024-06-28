from typing import List, Dict, Tuple
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
PARENT = 'Parent'
LABEL = 'Label'
COUNT = 'Count'
REF_COUNT = 'Reference count'
PROP = 'Proportion'
REF_PROP = 'Reference proportion'
RELAT_PROP = 'Relative proportion'
PVAL = 'Pvalue'


# ==================================================================================================
# FUNCTIONS
# ==================================================================================================

# Generate Data table
# --------------------------------------------------------------------------------------------------
def get_fig_parameters(classes_abondance: Dict[str, int], parent_dict: Dict[str, List[str]],
                       children_dict: Dict[str, List[str]], root_item,
                       subset_abundance: Dict[str, int] = None, full: bool = True,
                       names: Dict[str, str] = None) -> Dict[str, List]:
    """ Returns a dictionary of parameters to create the sunburst figure.

    Parameters
    ----------
    classes_abondance: Dict[str, int]
        Dictionary associating for each class the number of metabolites found belonging to the class.
    parent_dict: Dict[str, List[str]]
        Dictionary associating for each class, its parents classes
    children_dict: Dict[str, List[str]]
        Dictionary associating for each class, its children classes
    full: bool (default=True)
        True for a full figure with class duplication if a class has +1 parents.
        False to have a reduced vew with all label appearing only once (1 random parent chosen)
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
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
            - reference abundance value : Reference_count (int)
    """
    data = {IDS: list(),
            PARENT: list(),
            LABEL: list(),
            COUNT: list(),
            REF_COUNT: list()}
    for c_label, c_abundance in classes_abondance.items():
        c_sub_abundance = get_sub_abundance(subset_abundance, c_label, c_abundance)

        if c_label != root_item:
            c_parents = parent_dict[c_label]
            m_id = c_label
            label = c_label

            if names is not None:
                try:
                    label = names[c_label]
                except KeyError:
                    label = c_label

            data = add_value_data(data=data,
                                  m_id=m_id,
                                  label=label,
                                  value=c_sub_abundance,
                                  base_value=c_abundance,
                                  parent=c_parents[0])

            if len(c_parents) > 1 and full:
                for p in c_parents[1:]:
                    suffix = '__' + p
                    c_id = c_label + suffix
                    data = add_value_data(data=data,
                                          m_id=c_id,
                                          label=c_label,
                                          value=c_sub_abundance,
                                          base_value=c_abundance,
                                          parent=p)

                    if c_label in children_dict.keys():
                        c_children = children_dict[c_label]
                        for c in c_children:
                            data = add_children(data=data,
                                                origin=suffix,
                                                child=c,
                                                parent=c_id,
                                                classes_abondance=classes_abondance,
                                                subset_abundance=subset_abundance,
                                                children_dict=children_dict)
        else:
            data = add_value_data(data=data,
                                  m_id=c_label,
                                  label=c_label,
                                  value=c_sub_abundance,
                                  base_value=c_abundance,
                                  parent='')

    return data


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


def add_value_data(data: Dict[str, List], m_id: str, label: str, value: float, base_value: float,
                   parent: str) -> Dict[str, List]:
    """ Fill the data dictionary for a metabolite class.

    Parameters
    ----------
    data: Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (float)
            - reference abundance value : Reference_count (float)
    m_id: str
        ID of the metabolite class to add
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
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
            - reference abundance value : Reference_count (int)
    """
    if m_id in data[IDS]:
        raise ValueError(f'{m_id} already in data IDs, all IDs must be unique.')
    data[IDS].append(m_id)
    data[LABEL].append(label)
    data[PARENT].append(parent)
    data[COUNT].append(value)
    data[REF_COUNT].append(base_value)
    return data


def add_children(data: Dict[str, List], origin: str, child: str, parent: str,
                 classes_abondance: Dict[str, int],
                 subset_abundance, children_dict: Dict[str, List[str]]) \
        -> Dict[str, List]:
    """ Add recursively all children of a given class to the data dictionary.

    Parameters
    ----------
    data: Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
            - reference abundance value : Reference_count (int)
    origin: str
        Origin of propagation : parent class of parent
    child: str
        Child metabolite class
    parent: str
        Parent metabolite class
    classes_abondance: Dict[str, int]
        Dictionary associating for each class the number of metabolites found belonging to the class.
    subset_abundance
    children_dict: Dict[str, List[str]]
        Dictionary associating for each class, its children classes

    Returns
    -------
    Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
            - labels : Label (str)
            - parents ids : Parent (str)
            - abundance value : Count (int)
            - reference abundance value : Reference_count (int)
    """
    if child in classes_abondance.keys():
        c_sub_value = get_sub_abundance(subset_abundance, child, classes_abondance[child])

        data = add_value_data(data=data,
                              m_id=child + origin,
                              label=child,
                              value=c_sub_value,
                              base_value=classes_abondance[child],
                              parent=parent)

        if child in children_dict.keys():
            origin_2 = origin + '__' + child
            cs = children_dict[child]
            for c in cs:
                add_children(data=data,
                             origin=origin_2,
                             child=c,
                             parent=child + origin,
                             classes_abondance=classes_abondance,
                             subset_abundance=subset_abundance,
                             children_dict=children_dict)
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
        # TODO : IDK WHY IT WORKS ???
        missed = [data[IDS][i] for i in range(len(data[IDS])) if data[RELAT_PROP][i] < 1]
        if missed:
            parents = {data[PARENT][data[IDS].index(m)] for m in missed}
            for p in parents:
                data = get_relative_prop(data, p)
            missed = [data[IDS][i] for i in range(len(data[IDS])) if data[RELAT_PROP][i] < 1]
    return data


def get_relative_prop(data: Dict[str, List], p_id: str):
    """ Get recursively relative proportion of a parent children to itself. Add id to data
    Relative_proportion.

    Parameters
    ----------
    data: Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
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
                                 test: str, names: bool) \
        -> Tuple[Dict[str, List], Dict[str, float]]:
    """ Performs statistical tests for enrichment analysis.

    Parameters
    ----------
    data: Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
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
    names: bool
        True if names associated with labels, False otherwise

    Returns
    -------
    Dict[str, List]
        Dictionary with lists of :
            - ids : ID (str)
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
    M = np.max(list(ref_classes_abundance.values()))
    if names:
        m_list = [ref_classes_abundance[x] if x in ref_classes_abundance.keys() else 0 for x in
                  data[IDS]]
    else:
        m_list = [ref_classes_abundance[x] if x in ref_classes_abundance.keys() else 0 for x in
                  data[LABEL]]
    N = int(np.nanmax(data[COUNT]))
    n_list = data[COUNT]
    data[PVAL] = list()
    nb_classes = len(set([data[LABEL][i]
                          for i in range(len(data[COUNT]))
                          if data[COUNT][i] != np.nan]))
    significant_representation = dict()
    for i in range(len(m_list)):
        if type(n_list[i]) == int:
            # Binomial Test
            if test == BINOMIAL_TEST:
                p_val = stats.binomtest(n_list[i], N, m_list[i] / M, alternative='two-sided').pvalue
            # Hypergeometric Test
            elif test == HYPERGEO_TEST:
                p_val = stats.hypergeom.sf(n_list[i] - 1, M, m_list[i], N)
            else:
                raise ValueError(f'test parameter must be in : {[BINOMIAL_TEST, HYPERGEO_TEST]}')
            if ((n_list[i] / N) - (m_list[i] / M)) > 0:
                data[PVAL].append(-np.log10(p_val))
            else:
                data[PVAL].append(np.log10(p_val))
            if p_val < 0.05 / nb_classes:
                significant_representation[data[LABEL][i]] = p_val.round(10)
        else:
            data[PVAL].append(np.nan)
    significant_representation = dict(
        sorted(significant_representation.items(), key=lambda item: item[1]))
    return data, significant_representation