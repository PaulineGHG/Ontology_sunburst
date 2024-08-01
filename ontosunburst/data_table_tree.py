import copy
from typing import List, Dict, Set
import numpy as np
from numpy import nan
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

# Path cut
PATH_UNCUT = 'uncut'
PATH_DEEPER = 'deeper'
PATH_HIGHER = 'higher'
PATH_BOUND = 'bound'


# ==================================================================================================
# CLASS
# ==================================================================================================
class DataTable:
    def __init__(self):
        self.ids = list()
        self.onto_ids = list()
        self.labels = list()
        self.parents = list()
        self.count = list()
        self.ref_count = list()
        self.prop = list()
        self.ref_prop = list()
        self.relative_prop = list()
        self.p_val = list()
        self.len = 0

    def __str__(self):
        string = ''
        data = self.get_data_dict()
        for k, v in data.items():
            string += f'{k}\n{"-"*len(k)}\n{v}\n'
        return string

    def get_data_dict(self):
        return {IDS: self.ids, ONTO_ID: self.onto_ids, LABEL: self.labels, PARENT: self.parents,
                COUNT: self.count, REF_COUNT: self.ref_count, PROP: self.prop,
                REF_PROP: self.ref_prop, RELAT_PROP: self.relative_prop, PVAL: self.p_val}

    def fill_parameters(self, ref_abundance: Dict[str, int], parent_dict: Dict[str, List[str]],
                        root_item: str, set_abundance: Dict[str, int],
                        names: Dict[str, str] = None, ref_base: bool = True):
        """ Fill DataTable list attributes (self.ids, self.onto_ids, self.labels, self.parents,
        self.count, self.ref_count)

        Parameters
        ----------
        ref_abundance: Dict[str, int]
            Dictionary associating for each class the number of objects found belonging to the class
        parent_dict: Dict[str, List[str]]
            Dictionary associating for each class, its parents classes
        root_item: str
            Name of the root item of the ontology
        set_abundance: Dict[str, int]
        names: Dict[str, str]
            Dictionary associating metabolic object ID to its Name
        ref_base
        """
        if ref_base:
            for c_onto_id, c_ref_abundance in ref_abundance.items():
                c_abundance = get_sub_abundance(set_abundance, c_onto_id, c_ref_abundance)
                self.__fill_id_parameter(c_onto_id, root_item, names, parent_dict, c_abundance,
                                         c_ref_abundance)
        else:
            for c_onto_id, c_abundance in set_abundance.items():
                c_ref_abundance = ref_abundance[c_onto_id]
                self.__fill_id_parameter(c_onto_id, root_item, names, parent_dict, c_abundance,
                                         c_ref_abundance)

    def __fill_id_parameter(self, c_onto_id: str, root_item: str, names: Dict[str, str],
                            parent_dict: Dict[str, List[str]], c_abundance: float,
                            c_ref_abundance: float):
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
                self.add_value(m_id=c_id, onto_id=c_onto_id, label=c_label,
                               count=c_abundance, ref_count=c_ref_abundance,
                               parent=c_id[len(c_onto_id) + 2:])  # Remove c_label__ prefix
        else:
            self.add_value(m_id=c_onto_id, onto_id=c_onto_id, label=c_onto_id,
                           count=c_abundance, ref_count=c_ref_abundance, parent='')

    def add_value(self, m_id: str, onto_id: str, label: str, count: float, ref_count: float,
                  parent: str):
        """ Fill the data attributes for an object class.

        Parameters
        ----------
        m_id: str
            ID unique of the object class to add
        onto_id: str
            ID in the ontology
        label: str
            Label (name) of the object class to add
        count: int
            Abundance value of the object class to add
        ref_count: int
            Reference abundance value of the object class to add
        parent: str
            Parent object class of the object class to add
        """
        if m_id in self.ids:
            raise ValueError(f'{m_id} already in data IDs, all IDs must be unique.')
        self.ids.append(m_id)
        self.onto_ids.append(onto_id)
        self.labels.append(label)
        self.parents.append(parent)
        self.count.append(count)
        self.ref_count.append(ref_count)
        self.prop.append(nan)
        self.ref_prop.append(nan)
        self.relative_prop.append(nan)
        self.p_val.append(nan)
        self.len += 1

    def calculate_proportions(self, ref_base: bool):
        """ Calculate DataTable proportion list attributes (self.prop, self.ref_prop,
        self.relative_prop). If total add relative proportion to +1 parent for branch value.

        Parameters
        ----------
        ref_base: bool
            True if reference base representation
        """
        # Get total proportion
        max_abondance = int(np.nanmax(self.count))
        self.prop = [x / max_abondance for x in self.count]
        # Get reference proportion
        max_ref_abondance = np.max(self.ref_count)
        self.ref_prop = [x / max_ref_abondance for x in self.ref_count]
        # Get proportion relative to +1 parent proportion for total branch value
        self.relative_prop = [x for x in self.prop]
        p = ''
        self.__get_relative_prop(p, ref_base)
        # IDK WHY IT WORKS ???
        missed = [self.ids[i] for i in range(self.len) if self.relative_prop[i] < 1]
        if missed:
            parents = {self.parents[self.ids.index(m)] for m in missed}
            for p in parents:
                self.__get_relative_prop(p, ref_base)

    def __get_relative_prop(self, p_id: str, ref_base: bool):
        """ Get recursively relative proportion of a parent children to itself. Set it to class
        self.relative_prop attribute.

        Parameters
        ----------
        p_id: str
            ID of the parent
        """
        if ref_base:
            base_count = self.ref_count
        else:
            base_count = self.count
        if p_id == '':
            prop_p = MAX_RELATIVE_NB
            count_p = max(base_count)
        else:
            prop_p = self.relative_prop[self.ids.index(p_id)]
            count_p = base_count[self.ids.index(p_id)]
        index_p = [i for i, v in enumerate(self.parents) if v == p_id]
        p_children = [self.ids[i] for i in index_p]
        count_p_children = [base_count[i] for i in index_p]
        if np.nansum(count_p_children) > count_p:
            total = np.nansum(count_p_children)
        else:
            total = count_p
        for i, c in enumerate(p_children):
            if not ref_base and np.isnan(self.prop[self.ids.index(c)]):
                prop_c = 0
            else:
                prop_c = int((count_p_children[i] / total) * prop_p)
            self.relative_prop[self.ids.index(c)] = prop_c
        for c in p_children:
            if c in self.parents:
                self.__get_relative_prop(c, ref_base)

    def make_enrichment_analysis(self, test: str, scores=None) -> Dict[str, float]:
        """ Performs statistical tests for enrichment analysis.

        Parameters
        ----------
        test: str
            Type of test : binomial or hypergeometric
        scores: Dict

        Returns
        -------
        Dict[str, float]
            Dictionary of significant metabolic object label associated with their p-value
        """
        nb_classes = len(set([self.labels[i] for i in range(self.len)
                              if not np.isnan(self.count[i])]))
        significant_representation = dict()
        if scores is not None:
            for i in range(self.len):
                p_val = scores[self.onto_ids[i]]
                self.p_val[i] = -np.log10(p_val)
                if p_val < 0.05 / nb_classes:  # Keep significant p-values : Bonferroni
                    significant_representation[self.onto_ids[i]] = p_val
        else:
            m = np.max(self.ref_count)  # M = ref set total item number
            n = int(np.nanmax(self.count))  # N = interest set total item number
            for i in range(self.len):
                if type(self.count[i]) == int:  # If count not nan (= if concept in interest set)
                    # Binomial Test
                    if test == BINOMIAL_TEST:
                        p_val = stats.binomtest(self.count[i], n, self.ref_count[i] / m,
                                                alternative='two-sided').pvalue
                    # Hypergeometric Test
                    elif test == HYPERGEO_TEST:
                        p_val_upper = stats.hypergeom.sf(self.count[i] - 1, m, self.ref_count[i], n)
                        p_val_lower = stats.hypergeom.cdf(self.count[i], m, self.ref_count[i], n)
                        p_val = 2 * min(p_val_lower, p_val_upper)  # bilateral
                    else:
                        raise ValueError(
                            f'test parameter must be in : {[BINOMIAL_TEST, HYPERGEO_TEST]}')
                    if ((self.count[i] / n) - (self.ref_count[i] / m)) > 0:  # If over-represented :
                        self.p_val[i] = -np.log10(p_val)  # Positive log10(p-value)
                    else:  # If under-represented :
                        self.p_val[i] = np.log10(p_val)  # Negative log10(p-value)
                    if p_val < 0.05 / nb_classes:  # Keep significant p-values : Bonferroni
                        significant_representation[self.onto_ids[i]] = p_val
        significant_representation = dict(
            sorted(significant_representation.items(), key=lambda item: item[1]))
        return significant_representation

    def cut_root(self, mode: str):
        """ Filter data to cut (or not) the root to remove not necessary 100% represented classes.

        Parameters
        ----------
        mode: str
            Mode of root cutting
            - uncut: doesn't cut and keep all nodes from ontology root
            - cut: keep only the lowest level 100% shared node
            - total: remove all 100% shared nodes (produces a pie at center)
        """
        if mode not in {ROOT_UNCUT, ROOT_CUT, ROOT_TOTAL_CUT}:
            raise ValueError(f'Root cutting mode {mode} unknown, '
                             f'must be in {[ROOT_UNCUT, ROOT_CUT, ROOT_TOTAL_CUT]}')
        if mode == ROOT_CUT or mode == ROOT_TOTAL_CUT:
            roots_ind = [i for i in range(self.len) if self.relative_prop[i] == MAX_RELATIVE_NB]
            roots = [self.ids[i] for i in roots_ind]
            self.delete_value(roots_ind)
            if mode == ROOT_CUT:
                self.parents = [str(x).split('__')[0] if x in roots else x for x in self.parents]
            if mode == ROOT_TOTAL_CUT:
                self.parents = ['' if x in roots else x for x in self.parents]

    def cut_nested_path(self, mode: str):
        if mode != PATH_UNCUT:
            nested_paths = []
            for p_i in range(self.len):
                p = self.ids[p_i]
                p_children = [self.ids[i] for i in range(self.len) if self.parents[i] == p]
                if len(p_children) == 1:
                    p_p = self.parents[p_i]
                    p_p_children = [self.ids[i] for i in range(self.len) if self.parents[i] == p_p]
                    if len(p_p_children) != 1:
                        p_count = self.count[p_i]
                        c_i = self.ids.index(p_children[0])
                        c_count = self.count[c_i]
                        if p_count == c_count:
                            nested_paths.append(self.get_full_nested_path(c_i, [p_i]))
            self.delete_nested_path(mode, nested_paths)

    def get_full_nested_path(self, p_i, n_path):
        n_path.append(p_i)
        p = self.ids[p_i]
        p_children = [self.ids[i] for i in range(self.len) if self.parents[i] == p]
        if len(p_children) == 1:
            p_count = self.count[p_i]
            c_i = self.ids.index(p_children[0])
            c_count = self.count[c_i]
            if p_count == c_count:
                n_path = self.get_full_nested_path(c_i, n_path)
        return n_path

    def delete_nested_path(self, mode, nested_paths):
        to_del = []
        if mode == PATH_DEEPER:
            for path in nested_paths:
                to_del += path[:-1]
                to_keep = path[-1]
                root_p = self.parents[path[0]]
                self.parents[to_keep] = root_p
                self.labels[to_keep] = '... ' + self.labels[to_keep]
        elif mode == PATH_HIGHER:
            for path in nested_paths:
                to_del += path[1:]
                to_keep = path[0]
                to_keep_c = [i for i in range(self.len) if self.parents[i] == self.ids[path[-1]]]
                for c_i in to_keep_c:
                    self.parents[c_i] = self.ids[to_keep]
                self.labels[to_keep] += ' ...'
        elif mode == PATH_BOUND:
            for path in nested_paths:
                to_del += path[1:-1]
                to_keep_up = path[0]
                to_keep_do = path[-1]
                self.parents[to_keep_do] = self.ids[to_keep_up]
                self.labels[to_keep_up] += ' ...'
                self.labels[to_keep_do] = '... ' + self.labels[to_keep_do]
        self.delete_value(to_del)

    def delete_value(self, v_index: int or List[int]):
        data = self.get_data_dict()
        if type(v_index) == int:
            v_index = [v_index]
        for i in sorted(v_index, reverse=True):
            for k, v in data.items():
                del v[i]
            self.len -= 1

    def get_col(self, index: int or List[int] = None) -> List or List[List]:
        if index is None:
            index = list(range(self.len))
        if type(index) == int:
            index = [index]
        cols = list()
        for i in index:
            cols.append((self.ids[i], self.onto_ids[i], self.labels[i], self.parents[i],
                         self.count[i], self.ref_count[i], self.prop[i], self.ref_prop[i],
                         self.relative_prop[i], self.p_val[i]))
        return cols


# ==================================================================================================
# FUNCTIONS
# ==================================================================================================
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
