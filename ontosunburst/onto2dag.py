import logging
from typing import List, Set, Dict, Any, TypeAlias, Literal
import numpy
from ontosunburst.input_preprocessing import OntologyDAG, IdToLabel, InputsAb, Weight


class SubDAG:
    def __init__(self, interest: InputsAb, reference: InputsAb, all_concepts: set[str],
                 ontology_dag: OntologyDAG, root: str, id_to_labels: IdToLabel):
        self.nodes = []
        self.root = root
        self.leaves = []
        ontology_children_dag = get_children_dict(ontology_dag)
        concepts_ancestors = get_ancestors(all_concepts, ontology_dag, root)
        i_cum_w = get_cumulative_w(concepts_ancestors, interest)
        r_cum_w = get_cumulative_w(concepts_ancestors, reference)
        for c in all_concepts:
            node = NodeDAG(onto_id=c,
                           label=dict_value_or(id_to_labels, c, c),
                           exp_w=dict_value_or(interest, c, numpy.nan),
                           cum_w=dict_value_or(i_cum_w, c, numpy.nan),
                           max_w=i_cum_w[root],
                           r_exp_w=dict_value_or(reference, c, numpy.nan),
                           r_cum_w=dict_value_or(r_cum_w, c, numpy.nan),
                           r_max_w=r_cum_w[root],
                           parents=dict_value_or(ontology_dag, c, []),
                           children=dict_value_or(ontology_children_dag, c, []))
            self.nodes.append(node)


class NodeDAG:
    def __init__(self, onto_id: str, label: str, exp_w: Weight, cum_w: Weight, max_w: Weight,
                 r_exp_w: Weight, r_cum_w: Weight, r_max_w: Weight,
                 parents: List[str], children: List[str]):
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
        self.enrichment_p_val = None
        self.difference = cum_w - r_cum_w
        # Hierarchy
        self.parents = parents
        self.children = children


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
def get_cumulative_w(concepts_ancestors: Dict[str, Set[str]], inputs_ab: InputsAb) -> \
        Dict[str, float]:
    """ Calculate the cumulative weight of each class depending on the inputs abundances
    (interest or reference).

    Parameters
    ----------
    concepts_ancestors: Dict[str, Set[str]] (Dict[metabolite, Set[class]])
        Dictionary associating for each concept the list of all parent classes it belongs to.
    inputs_ab: Dict[str, float]
        Dictionary associating for each concept, its abundance value.

    Returns
    -------
    Dict[str, float]
        Dictionary associating for each class ist cumulative weight.
    """
    cumulative_weights = dict()
    for cpt, ab in inputs_ab.items():
        if cpt not in cumulative_weights.keys():
            cumulative_weights[cpt] = ab
        else:
            cumulative_weights[cpt] += ab
        for c in concepts_ancestors[cpt]:
            if c not in cumulative_weights.keys():
                cumulative_weights[c] = ab
            else:
                cumulative_weights[c] += ab
    return dict(reversed(sorted(cumulative_weights.items(), key=lambda item: item[1])))


def dict_value_or(dictionary: dict, value: str, alternative: Any):
    try:
        return dictionary[value]
    except KeyError:
        return alternative


def get_children_dict(parent_dict: OntologyDAG) -> OntologyDAG:
    """ Create the children dictionary from the parents dictionary.
    Parameters
    ----------
    parent_dict: dict[str, list[str]]
        Dictionary associating for each class, its parents classes
    Returns
    -------
    dict[str, list[str]]
        Dictionary associating for each class, its children classes
    """
    children_dict = dict()
    for c, ps in parent_dict.items():
        for p in ps:
            if p not in children_dict.keys():
                children_dict[p] = list()
            if c not in children_dict.keys():
                children_dict[c] = list()
            children_dict[p].append(c)
    return children_dict


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
