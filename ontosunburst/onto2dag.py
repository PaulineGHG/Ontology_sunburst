import logging
from typing import List, Set, Dict, Any, TypeAlias, Literal
import numpy
from ontosunburst.input_preprocessing import OntologyDAG, IdToLabel, InputsAb


class SubDAG:
    def __init__(self, interest: InputsAb, reference: InputsAb, all_concepts: set[str],
                 ontology_dag: OntologyDAG, root: str, id_to_labels: IdToLabel):
        self.nodes = []
        self.root = root
        self.leaves = []

        concepts_all_parents = get_all_parents(all_concepts, ontology_dag, root)

        # abundances_dict = get_abundance_dict(abundances, concepts)
        # ref_abundances_dict = get_abundance_dict(ref_abundances, ref_concepts)
        #
        # cum_w = get_cumulative_w(concepts_all_classes, abundances_dict)
        # r_cum_w = get_cumulative_w(concepts_all_classes, ref_abundances_dict)
        #
        # for c in concepts:
        #     if c in ref_abundances_dict:
        #         node = NodeDAG(onto_id=c, label=id_to_labels[c],
        #                        exp_w=abundances_dict[c], cum_w=cum_w[c], max_w=cum_w[root],
        #                        r_exp_w=ref_abundances_dict[c], r_cum_w=r_cum_w[c],
        #                        r_max_w=r_cum_w[root])


class NodeDAG:
    def __init__(self, onto_id, label, exp_w, cum_w, max_w, r_exp_w, r_cum_w, r_max_w):
        # ID and label
        self.onto_id = onto_id
        self.label = label
        # Weights
        self.experimental_weight = exp_w
        self.cumulative_weight = cum_w
        self.proportion = cum_w / max_w
        # Reference weights
        self.ref_experimental_weight = r_exp_w
        self.ref_cumulative_weight = r_cum_w
        self.ref_proportion = r_cum_w / r_max_w
        # Comparison calculations
        self.intensity = None
        self.difference = cum_w - r_cum_w
        # Hierarchy
        self.parents = []
        self.children = []


# Main ontology to reduced dag functions
# --------------------------------------------------------------------------------------------------
def reduce_d_ontology(complete_dictionary: Dict[str, Any],
                      classes_abundance: Dict[str, float]) -> Dict[str, Any]:
    """ Extract the sub-graph of the d_classes_ontology dictionary conserving only nodes implicated
    with the concepts studied.

    Parameters
    ----------
    complete_dictionary: Dict[str, Any]
        Dictionary of the ontology complete graph
    classes_abundance: Dict[str, float]
        Dictionary of abundances (keys are all nodes implicated to be conserved)

    Returns
    -------
    Dict[str, Any]
        Dictionary of the ontology sub-graph conserving only nodes implicated with the concepts
        studied.
    """
    if complete_dictionary is not None:
        reduced_dictionary = dict()
        for k, v in complete_dictionary.items():
            if k in classes_abundance:
                reduced_dictionary[k] = v
        return reduced_dictionary


# ==================================================================================================
# REDUCE DAG FUNCTIONS
# ==================================================================================================


# Recursive class extraction function
# --------------------------------------------------------------------------------------------------
def get_all_parents(all_concepts: Set[str], ontology_dag: OntologyDAG, root: str) \
        -> Dict[str, Set[str]]:
    all_parents_dict = dict()
    # for met, classes in obj_classes.items():
    #     all_classes = set(classes)
    #     for c in classes:
    #         if c != root_item:
    #             m_classes = get_parents(c, set(d_classes_ontology[c]), d_classes_ontology,
    #                                     root_item)
    #             all_classes = all_classes.union(m_classes)
    #     all_parents[met] = all_classes
    for cpt in all_concepts:
        parents = set(ontology_dag[cpt])
        all_parents = get_parents(cpt, parents, ontology_dag, root)
        all_parents_dict[cpt] = all_parents
    return all_parents_dict


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
# WEIGHTS CALCULATION
# ==================================================================================================

def get_abundance_dict(abundances: List[float] or None, concepts: List[str]) \
        -> Dict[str, float]:
    """ Generate abundances dictionary.

    Parameters
    ----------
    abundances: List[float] (size N) or None
        List of concepts abundances (or None if no abundances associated --> will associate
        an abundance of 1 for each concept)
    concepts: List[str] (size N)
        List of concepts ID.

    Returns
    -------
    Dict[str, float]
        Dictionary associating to each concept its abundance.
    """
    if abundances is None:
        abundances = len(concepts) * [1]
    if len(concepts) == len(abundances):
        abundances_dict = {}
        for i in range(len(concepts)):
            abundances_dict[concepts[i]] = abundances[i]
    else:
        raise AttributeError(f'Length of concepts IDs list must be equal to '
                             f'its abundances list length : {len(concepts)} '
                             f'!= {len(abundances)}')
    return abundances_dict


def calculate_weights(all_classes: Dict[str, Set[str]], abundances_dict: Dict[str, float],
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


def get_cumulative_w(all_classes: Dict[str, Set[str]], abundances_dict: Dict[str, float]) -> Dict[str, float]:
    """ Indicate for each class the number of base object found belonging to the class

    Parameters
    ----------
    all_classes: Dict[str, Set[str]] (Dict[metabolite, Set[class]])
        Dictionary associating for each concept the list of all parent classes it belongs to.
    abundances_dict: Dict[str, float]
        Dictionary associating for each concept, its abundance value

    Returns
    -------
    Dict[str, float]
        Dictionary associating for each class the weight of concepts found belonging to the class.
    """
    classes_abondance = dict()
    for met, ab in abundances_dict.items():
        if met not in classes_abondance.keys():
            classes_abondance[met] = ab
        else:
            classes_abondance[met] += ab
        for c in all_classes[met]:
            if c not in classes_abondance.keys():
                classes_abondance[c] = ab
            else:
                classes_abondance[c] += ab
    return dict(reversed(sorted(classes_abondance.items(), key=lambda item: item[1])))


def get_classes_scores(all_classes, scores_dict, root):
    if scores_dict is not None:
        classes_scores = dict()
        for met, classes in all_classes.items():
            if met in scores_dict.keys():
                classes_scores[met] = scores_dict[met]
            else:
                classes_scores[met] = numpy.nan
        classes_scores[root] = numpy.nan
        return classes_scores
