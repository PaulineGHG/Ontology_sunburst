from typing import List, Set, Dict, Any
import numpy


# ==================================================================================================
# CLASSES EXTRACTION
# ==================================================================================================

# Main class extraction function
# --------------------------------------------------------------------------------------------------
def extract_classes(concepts: List[str], root: str,
                    d_classes_ontology: Dict[str, List[str]]) -> Dict[str, Set[str]]:
    """ Extract all parent classes (until root) from a list of concepts.

    Parameters
    ----------
    concepts: List[str]
        List of concepts to classify
    root: str
        Root item of the ontology
    d_classes_ontology: Dict[str, List[str]]
        Dictionary of the classes ontology associating for each concept its +1 parent classes.

    Returns
    -------
    Dict[str, Set[str]]
        Dictionary associating to each concept all its parent classes (until the root)
    """
    d_obj_classes = dict()
    print(f'{len(concepts)} concepts to classify')
    for obj in concepts:
        try:
            d_obj_classes[obj] = d_classes_ontology[obj]
            classified = True
        except KeyError:
            classified = False
        if not classified:
            print(f'{obj} not classified.')
    print(f'{len(d_obj_classes)}/{len(concepts)} concepts classified')
    return get_all_classes(d_obj_classes, d_classes_ontology, root)


# Recursive class extraction function
# --------------------------------------------------------------------------------------------------
def get_all_classes(obj_classes: Dict[str, List[str]], d_classes_ontology: Dict[str, List[str]],
                    root_item: str) -> Dict[str, Set[str]]:
    """ Extract all parent classes for each metabolite.

    Parameters
    ----------
    obj_classes: Dict[str, List[str]] (Dict[metabolite, List[class]])
        Dictionary associating for each object the list of +1 parent classes it belongs to.
    d_classes_ontology: Dict[str, List[str]]
        Dictionary of the classes ontology associating for each class its +1 parent classes.
    root_item: str
        Name of the root item of the ontology.

    Returns
    -------
    Dict[str, Set[str]] (Dict[metabolite, Set[class]])
        Dictionary associating for each metabolite the list of all parent classes it belongs to.
    """
    all_classes_met = dict()
    for met, classes in obj_classes.items():
        all_classes = set(classes)
        for c in classes:
            if c != root_item:
                m_classes = get_parents(c, set(d_classes_ontology[c]), d_classes_ontology, root_item)
                all_classes = all_classes.union(m_classes)
        all_classes_met[met] = all_classes
    return all_classes_met


def get_parents(child: str, parent_set: Set[str], d_classes_ontology: Dict[str, List[str]],
                root_item) -> Set[str]:
    """ Get recursively from a child class, all its parents classes found in ontology.

    Parameters
    ----------
    child: str
        Child class
    parent_set: Set[str]
        Set of all parents from previous classes
    d_classes_ontology: Dict[str, List[str]]
        Dictionary of the classes ontology of MetaCyc associating for each class its parent classes.
    root_item: str
        Name of the root item of the ontology

    Returns
    -------
    Set[str]
        Set of the union of the set  of child parent classes and the set of all previous parents.
    """
    parents = d_classes_ontology[child]
    for p in parents:
        parent_set.add(p)
    for p in parents:
        if p != root_item:
            parent_set = get_parents(p, parent_set, d_classes_ontology, root_item)
    return parent_set


# ==================================================================================================
# ABUNDANCES CALCULATION
# ==================================================================================================

def get_abundance_dict(abundances: List[float] or None, metabolic_objects: List[str], ref: bool)\
        -> Dict[str, float]:
    """ Generate abundances dictionary.

    Parameters
    ----------
    abundances: List[float] (size N) or None
        List of metabolic objects abundances (or None if no abundances associated --> will associate
        an abundance of 1 for each object)
    metabolic_objects: List[str] (size N)
        List of metabolic objects ID.
    ref: bool
        True if metabolic objects are the reference list, False otherwise (subset / study case).

    Returns
    -------
    Dict[str, float]
        Dictionary associating to each object its abundance.
    """
    if abundances is None:
        abundances = len(metabolic_objects) * [1]
    if len(metabolic_objects) == len(abundances):
        abundances_dict = {}
        for i in range(len(metabolic_objects)):
            abundances_dict[metabolic_objects[i]] = abundances[i]
    else:
        if ref:
            raise AttributeError(f'Length of "reference_set" parameter must be equal to '
                                 f'"ref_abundances" parameter length : {len(metabolic_objects)} '
                                 f'!= {len(abundances)}')
        else:
            raise AttributeError(f'Length of "metabolic_objects" parameter must be equal to '
                                 f'"abundances" parameter length : {len(metabolic_objects)} '
                                 f'!= {len(abundances)}')
    return abundances_dict


def get_classes_abundance(all_classes: Dict[str, Set[str]], abundances_dict: Dict[str, float],
                          show_leaves: bool) -> Dict[str, float]:
    """ Indicate for each class the number of base object found belonging to the class

    Parameters
    ----------
    all_classes: Dict[str, Set[str]] (Dict[metabolite, Set[class]])
        Dictionary associating for each concept the list of all parent classes it belongs to.
    abundances_dict: Dict[str, float]
        Dictionary associating for each concept, its abundance value
    show_leaves: bool
        True to show input metabolic objets at sunburst leaves

    Returns
    -------
    Dict[str, float]
        Dictionary associating for each class the weight of concepts found belonging to the class.
    """
    classes_abondance = dict()
    for met, classes in all_classes.items():
        if show_leaves:
            if met not in classes_abondance.keys():
                classes_abondance[met] = abundances_dict[met]
            else:
                classes_abondance[met] += abundances_dict[met]
        for c in classes:
            if c not in classes_abondance.keys():
                classes_abondance[c] = abundances_dict[met]
            else:
                classes_abondance[c] += abundances_dict[met]
    return dict(reversed(sorted(classes_abondance.items(), key=lambda item: item[1])))


def get_classes_scores(all_classes, scores_dict, root):
    classes_scores = dict()
    for met, classes in all_classes.items():
        if met in scores_dict.keys():
            classes_scores[met] = scores_dict[met]
        else:
            classes_scores[met] = numpy.nan
    classes_scores[root] = numpy.nan
    return classes_scores


# ==================================================================================================
# UTILS
# ==================================================================================================

def reduce_d_ontology(d_classes_ontology: Dict[str, Any],
                      classes_abundance: Dict[str, float]) -> Dict[str, Any]:
    """ Extract the sub-graph of the d_classes_ontology dictionary conserving only nodes implicated
    with the concepts studied.

    Parameters
    ----------
    d_classes_ontology: Dict[str, Any]
        Dictionary of the ontology complete graph
    classes_abundance: Dict[str, float]
        Dictionary of abundances (keys are all nodes implicated to be conserved)

    Returns
    -------
    Dict[str, Any]
        Dictionary of the ontology sub-graph conserving only nodes implicated with the concepts
        studied.
    """
    reduced_d_ontology = dict()
    for k, v in d_classes_ontology.items():
        if k in classes_abundance:
            reduced_d_ontology[k] = v
    return reduced_d_ontology

