from typing import List, Tuple
from ontosunburst.onto2dag import NodeDAG, SubDAG

# ==================================================================================================
# CONSTANTS
# ==================================================================================================
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
class TreeNode:

    def __init__(self, dag_node: NodeDAG, parent, copy: int, r_prop: float):
        self.dag_node: NodeDAG = dag_node
        self.parent: TreeNode | None = parent
        self.node_id: Tuple[str, int] = (self.dag_node.onto_id, copy)
        self.relative_proportion: float = r_prop

    def get_arguments(self):
        dag_node_arg = self.dag_node.get_arguments()
        parent = None
        if self.parent is not None:
            parent = self.parent.node_id
        tree_node_arg = {'ID': self.node_id,
                         'Parent': parent,
                         'Relative proportion': self.relative_proportion}
        return tree_node_arg


class InducedTree:

    def __init__(self, sub_dag: SubDAG, ref_base: bool):
        self.nodes = {}
        self.sub_dag = sub_dag

        root_id = sub_dag.root
        root_dag_node = sub_dag.node_by_id(root_id)
        r_children = root_dag_node.children
        root_tree_node = TreeNode(root_dag_node, None, 0, 1.0)
        self.nodes[(root_id, 0)] = root_tree_node
        self.dag_traversal_rec(root_tree_node, r_children, ref_base)

    def dag_traversal_rec(self, parent: TreeNode, children: List[str], ref_base: bool):
        p_r_prop = parent.relative_proportion
        child_nodes = self.sub_dag.nodes_by_id(children)
        c_nodes_prop_sum = sum([c.proportion for c in child_nodes])
        if ref_base:
            total = parent.dag_node.ref_proportion
        else:
            total = parent.dag_node.proportion
        if c_nodes_prop_sum > total:
            total = c_nodes_prop_sum
        for c in child_nodes:
            if ref_base:
                c_r_prop = (c.ref_proportion / total) * p_r_prop
            else:
                c_r_prop = (c.proportion / total) * p_r_prop
            c_tree_node = TreeNode(c, parent, c.copies, truncate(c_r_prop, 5))
            self.nodes[(c.onto_id, c.copies)] = c_tree_node
            c.copies += 1
            c_children = c.children
            if c_children:
                self.dag_traversal_rec(parent=c_tree_node, children=c_children, ref_base=ref_base)


def truncate(n: float, dec: int) -> float:
    return float(int(n * (10 ** dec)) / (10 ** dec))


#     def cut_root(self, mode: str):
#         """ Filter data to cut (or not) the root to remove not necessary 100% represented classes.
#
#         Parameters
#         ----------
#         mode: str
#             Mode of root cutting
#             - uncut: doesn't cut and keep all nodes from ontology root
#             - cut: keep only the lowest level 100% shared node
#             - total: remove all 100% shared nodes (produces a pie at center)
#         """
#         if mode not in {ROOT_UNCUT, ROOT_CUT, ROOT_TOTAL_CUT}:
#             raise ValueError(f'Root cutting mode {mode} unknown, '
#                              f'must be in {[ROOT_UNCUT, ROOT_CUT, ROOT_TOTAL_CUT]}')
#         if mode == ROOT_CUT or mode == ROOT_TOTAL_CUT:
#             roots_ind = [i for i in range(self.len) if self.relative_prop[i] == MAX_RELATIVE_NB]
#             roots = [self.ids[i] for i in roots_ind]
#             roots_lab = [self.labels[i] if self.labels[i] not in self.ids
#                          else self.labels[i] + '_' for i in roots_ind]
#             lab = {roots[i]: roots_lab[i] for i in range(len(roots))}
#             self.delete_value(roots_ind)
#             if mode == ROOT_CUT:
#                 self.parents = [lab[x] if x in roots else x for x in self.parents]
#             if mode == ROOT_TOTAL_CUT:
#                 self.parents = ['' if x in roots else x for x in self.parents]
#
#     def cut_nested_path(self, mode: str, ref_base: bool):
#         """ Cut nested path in the tree graph (path of nested sectors sharing the same value)
#
#         Parameters
#         ----------
#         mode: str
#             Mode of path cutting
#             - uncut: doesn't cut and keep all sectors
#             - deeper: cut nested path and only conserve the deepest sector in the tree
#             - higher: cut nested path and only conserve the highest sector in the tree
#             - bound: cut nested path and only conserve the highest AND the deepest sectors in the
#             tree
#         ref_base: bool
#             True if reference base representation
#         """
#         if ref_base:
#             count = self.ref_count
#         else:
#             count = self.count
#         if mode != PATH_UNCUT:
#             nested_paths = []
#             for p_i in range(self.len):
#                 p = self.ids[p_i]
#                 p_children = [self.ids[i] for i in range(self.len) if self.parents[i] == p]
#                 if len(p_children) == 1:
#                     p_p = self.parents[p_i]
#                     p_p_children = [self.ids[i] for i in range(self.len) if self.parents[i] == p_p]
#                     if len(p_p_children) != 1:
#                         p_count = count[p_i]
#                         c_i = self.ids.index(p_children[0])
#                         c_count = count[c_i]
#                         if p_count == c_count:
#                             nested_paths.append(self.get_full_nested_path(c_i, [p_i], count))
#             self.delete_nested_path(mode, nested_paths)
#
#     def get_full_nested_path(self, p_i: int, n_path: List[int], count: List[float]):
#         """ Get all index of a nested path sector from its parent sector index.
#
#         Parameters
#         ----------
#         p_i: int
#             Parent sector index of the nested path
#         n_path: List[int]
#             List of sector indexes of the nested path
#         count: List[float]
#             List of all sectors count value.
#
#         Returns
#         -------
#         List[int]
#             List of sector indexes of the nested path
#         """
#         n_path.append(p_i)
#         p = self.ids[p_i]
#         p_children = [self.ids[i] for i in range(self.len) if self.parents[i] == p]
#         if len(p_children) == 1:
#             p_count = count[p_i]
#             c_i = self.ids.index(p_children[0])
#             c_count = count[c_i]
#             if p_count == c_count:
#                 n_path = self.get_full_nested_path(c_i, n_path, count)
#         return n_path
#
#     def delete_nested_path(self, mode: str, nested_paths: List[List[int]]):
#         """ Delete some sectors of the nested path to conserve only the deepest (deeper mode), only
#         the highest (higher mode) or both (bound mode)
#
#         Parameters
#         ----------
#         mode: str
#             Mode of path cutting
#             - uncut: doesn't cut and keep all sectors
#             - deeper: cut nested path and only conserve the deepest sector in the tree
#             - higher: cut nested path and only conserve the highest sector in the tree
#             - bound: cut nested path and only conserve the highest AND the deepest sectors in the
#             tree
#         nested_paths: List[List[int]]
#             List of lists of nested path sectors indexes
#         """
#         to_del = []
#         if mode == PATH_DEEPER:
#             for path in nested_paths:
#                 to_del += path[:-1]
#                 to_keep = path[-1]
#                 root_p = self.parents[path[0]]
#                 self.parents[to_keep] = root_p
#                 self.labels[to_keep] = '... ' + self.labels[to_keep]
#         elif mode == PATH_HIGHER:
#             for path in nested_paths:
#                 to_del += path[1:]
#                 to_keep = path[0]
#                 to_keep_c = [i for i in range(self.len) if self.parents[i] == self.ids[path[-1]]]
#                 for c_i in to_keep_c:
#                     self.parents[c_i] = self.ids[to_keep]
#                 self.labels[to_keep] += ' ...'
#         elif mode == PATH_BOUND:
#             for path in nested_paths:
#                 to_del += path[1:-1]
#                 to_keep_up = path[0]
#                 to_keep_do = path[-1]
#                 self.parents[to_keep_do] = self.ids[to_keep_up]
#                 if len(path) > 2:
#                     self.labels[to_keep_up] += ' ...'
#                     self.labels[to_keep_do] = '... ' + self.labels[to_keep_do]
#         self.delete_value(to_del)
#
#     def delete_value(self, v_index: int or List[int]):
#         """ Delete a sector of TreeData from its index or a list of sectors from a list of indexes
#
#         Parameters
#         ----------
#         v_index: int or List[int]
#             Index or list of indexes of the sectors to delete
#         """
#         data = self.get_data_dict()
#         if type(v_index) == int:
#             v_index = [v_index]
#         for i in sorted(v_index, reverse=True):
#             for k, v in data.items():
#                 del v[i]
#             self.len -= 1
#
#     def get_col(self, index: int or List[int] = None) -> List or List[List]:
#         """ Get a TreeData column from its index or a list of columns from a list of indexes.
#         Column = all values associated with a sector.
#
#         Parameters
#         ----------
#         index: int or List[int]
#             Index or list of indexes of the sectors to get the column
#
#         Returns
#         -------
#         List or List[List]
#             Column or list of columns obtained from indexes
#         """
#         if index is None:
#             index = list(range(self.len))
#         if type(index) == int:
#             index = [index]
#         cols = list()
#         for i in index:
#             cols.append((self.ids[i], self.onto_ids[i], self.labels[i], self.parents[i],
#                          self.count[i], self.ref_count[i], self.prop[i], self.ref_prop[i],
#                          self.relative_prop[i], self.p_val[i]))
#         return cols

