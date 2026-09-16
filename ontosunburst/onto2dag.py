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

        concepts_ancestors = get_ancestors(all_concepts, ontology_dag, root)
        i_cum_w = get_cumulative_w(concepts_ancestors, interest)
        r_cum_w = get_cumulative_w(concepts_ancestors, reference)

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


# Recursive ancestors extraction functions
# --------------------------------------------------------------------------------------------------
def get_ancestors(all_concepts: Set[str], ontology_dag: OntologyDAG, root: str) \
        -> Dict[str, Set[str]]:
    """ Return a dictionary associating to each input concept the list of all its ancestors
    (parents + parents of parents... recursively to the root)

    Parameters
    ----------
    all_concepts: set[str]
        Set of all concepts from interest and reference
    ontology_dag: dict[str, list[str]]
    root: str
        Root of the ontology DAG

    Returns
    -------
    dict[str, set[str]]
        Dictionary associating to each input concept the list of all its ancestors
        (parents + parents of parents... recursively to the root)
    """
    ancestors_dict = dict()
    for cpt in all_concepts:
        parents = set(ontology_dag[cpt])
        ancestors = get_ancestors_recursively(cpt, parents, ontology_dag, root)
        ancestors_dict[cpt] = ancestors
    return ancestors_dict


def get_ancestors_recursively(child: str, parents_set: Set[str], ontology_dag: OntologyDAG,
                              root_item) -> Set[str]:
    """ Get recursively from a child ID, all its ancestors IDs found in ontology.

    Parameters
    ----------
    child: str
        Child ID
    parents_set: set[str]
        Set of all parents from previous child ID
    ontology_dag: dict[str, list[str]]
    root_item: str
        Root of the ontology DAG

    Returns
    -------
    Set[str]
        Set of the union of the set of child parent classes and the set of all previous parents.
    """
    parents = ontology_dag[child]
    for p in parents:
        parents_set.add(p)
    for p in parents:
        if p != root_item:
            parents_set = get_ancestors_recursively(p, parents_set, ontology_dag, root_item)
    return parents_set


# ==================================================================================================
# WEIGHTS CALCULATION
# ==================================================================================================
# def calculate_weights(all_classes: Dict[str, Set[str]], abundances_dict: Dict[str, float],
#                       show_leaves: bool) -> Dict[str, float]:
#     """ Indicate for each class the number of base object found belonging to the class
#
#     Parameters
#     ----------
#     all_classes: Dict[str, Set[str]] (Dict[metabolite, Set[class]])
#         Dictionary associating for each concept the list of all parent classes it belongs to.
#     abundances_dict: Dict[str, float]
#         Dictionary associating for each concept, its abundance value
#     show_leaves: bool
#         True to show input metabolic objets at sunburst leaves
#
#     Returns
#     -------
#     Dict[str, float]
#         Dictionary associating for each class the weight of concepts found belonging to the class.
#     """
#     classes_abondance = dict()
#     for met, classes in all_classes.items():
#         if show_leaves:
#             if met not in classes_abondance.keys():
#                 classes_abondance[met] = abundances_dict[met]
#             else:
#                 classes_abondance[met] += abundances_dict[met]
#         for c in classes:
#             if c not in classes_abondance.keys():
#                 classes_abondance[c] = abundances_dict[met]
#             else:
#                 classes_abondance[c] += abundances_dict[met]
#     return dict(reversed(sorted(classes_abondance.items(), key=lambda item: item[1])))


def get_cumulative_w(concepts_ancestors: Dict[str, Set[str]], inputs_ab: InputsAb) -> \
        Dict[str, float]:
    """ Calculate the cumulative weight of each class depending on the inputs abundances
    (interest or reference).

    Parameters
    ----------
    concepts_ancestors: Dict[str, Set[str]] (Dict[metabolite, Set[class]])
        Dictionary associating for each concept the list of all parent classes it belongs to.
    inputs_ab: Dict[str, float]
        Dictionary associating for each concept, its abundance value

    Returns
    -------
    Dict[str, float]
        Dictionary associating for each class the weight of concepts found belonging to the class.
    """
    classes_abondance = dict()
    for cpt, ab in inputs_ab.items():
        if cpt not in classes_abondance.keys():
            classes_abondance[cpt] = ab
        else:
            classes_abondance[cpt] += ab
        for c in concepts_ancestors[cpt]:
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
