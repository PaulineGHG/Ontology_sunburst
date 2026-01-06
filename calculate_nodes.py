import networkx as nx
import matplotlib.pyplot as plt

import numpy as np

V = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
     'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']


G=nx.gnp_random_graph(25,0.05,directed=True)
roots = []
leaves = []
for n in G.nodes:
    root = True
    leaf = True
    for e in G.edges:
        if e[0] == n:
            root = False
        if e[1] == n:
            leaf = False
    if root:
        roots.append(n)
    if leaf:
        leaves.append(n)

ROOT = 'r'
G.add_node(ROOT)
for r in roots:
    G.add_edge(r, ROOT)

for n in G.nodes:
    print(list(nx.all_simple_paths(G, n, ROOT)))






