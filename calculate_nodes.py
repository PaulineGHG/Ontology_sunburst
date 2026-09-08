import networkx as nx
from ontosunburst import ontosunburst
import matplotlib.pyplot as plt

import numpy as np

V = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
     'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']
ROOT = 'r'

def gen_random_dag(main_root, nodes_nb=25):
    graph=nx.gnp_random_graph(nodes_nb,0.05,directed=True)
    root_nodes = []
    for n in graph.nodes:
        root_n = True
        for e in graph.edges:
            if e[0] == n:
                root_n = False
        if root_n:
            root_nodes.append(n)
    graph.add_node(main_root)
    for r in root_nodes:
        graph.add_edge(r, main_root)
    graph = remove_cycles(graph)
    graph = add_no_connected_nodes(graph)
    cycles = nx.recursive_simple_cycles(graph)
    if cycles:
        graph = gen_random_dag(main_root)
    leaf_nodes = []
    for n in graph.nodes:
        leaf_n = True
        for e in graph.edges:
            if e[1] == n:
                leaf_n = False
        if leaf_n:
            leaf_nodes.append(str(n))
    return graph, leaf_nodes


def remove_cycles(graph):
    cycles = nx.recursive_simple_cycles(graph)
    for c in cycles:
        if (c[-1], c[0]) in list(graph.edges):
            graph.remove_edge(c[-1], c[0])
    return graph

def add_no_connected_nodes(graph):
    not_connected = []
    paths = []
    for n in graph.nodes:
        paths_to_root = list(nx.all_simple_paths(graph, n, ROOT))
        if not paths_to_root:
            not_connected.append(n)
        else:
            paths += paths_to_root
    paths.sort(key=len)
    for n in not_connected:
        graph.add_edge(n, paths[1][0])
        del paths[1]
    return graph


def graph_to_dict(graph, main_root):
    graph_dict = {str(x):[] for x in graph.nodes if x != main_root}
    for e in graph.edges:
        graph_dict[str(e[0])].append(str(e[1]))
    return graph_dict


def count_all_paths_to_root(graph):
    c = 0
    for n in graph.nodes:
        if n != ROOT:
            paths = list(nx.all_simple_paths(graph, n, ROOT))
            c += len(paths)
    return c

def gen_figure(nb_dag_nodes):
    g, leaves = gen_random_dag(ROOT, nb_dag_nodes)
    g_dict = graph_to_dict(g, ROOT)
    nx.draw(g, with_labels=True)
    fig = ontosunburst(interest_set=leaves, ontology_dag_input=g_dict, input_root=ROOT, show_leaves=True)
    nb_sb_nodes = len(fig.to_dict()['data'][0]['ids'])
    nb_dag_paths = count_all_paths_to_root(g)
    print(nb_dag_nodes, nb_sb_nodes, nb_dag_paths)


gen_figure(20)


