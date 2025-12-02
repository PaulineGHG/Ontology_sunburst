from ontosunburst import ontosunburst
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def reverse_dag(original_dag):
    reversed_dag = {}
    for n, links in original_dag.items():
        for l in links:
            if l not in reversed_dag:
                reversed_dag[l] = []
            reversed_dag[l].append(n)
    return reversed_dag


def get_ontology_sub_dag(new_root, onto_dag):
    return reverse_dag(get_sub_dag(new_root, reverse_dag(onto_dag), {}))


def get_sub_dag(root, c_dag, sub_dag):
    sub_dag[root] = c_dag[root]
    for c in c_dag[root]:
        if c in c_dag:
            sub_dag = get_sub_dag(c, c_dag, sub_dag)
    return sub_dag


def load_onto_file(onto_file):
    with open(onto_file, 'r') as f:
        onto_dict = json.load(f)
    return onto_dict


def generate_ontology_sunburst(onto_dict, onto):
    lv_classes = set(onto_dict.keys())
    for p_classes in onto_dict.values():
        for p_class in p_classes:
            if p_class in lv_classes:
                lv_classes.remove(p_class)
    print(len(lv_classes))
    ontosunburst(interest_set=list(lv_classes), ontology=onto, show_leaves=True,
                 )


def decompose_metacyc(onto_dict):
    metacyc_reaction = get_ontology_sub_dag('Reactions', onto_dict)
    metacyc_pathway = get_ontology_sub_dag('Pathways', onto_dict)
    metacyc_compound = get_ontology_sub_dag('Compounds', onto_dict)
    metacyc_gene_product = get_ontology_sub_dag('Macromolecules', onto_dict)
    return metacyc_reaction, metacyc_pathway, metacyc_compound, metacyc_gene_product


def decompose_kegg(onto_dict):
    # PATHWAY
    kegg_modules = get_ontology_sub_dag('br:ko00002', onto_dict)
    kegg_pathway = get_ontology_sub_dag('br:br08901', onto_dict)
    kegg_pathway.update(kegg_modules)
    # REACTION
    kegg_reaction = get_ontology_sub_dag('br:br08201', onto_dict)
    # GENE PRODUCT
    kegg_gene_product = get_ontology_sub_dag('br:ko00001', onto_dict)
    # COMPOUND
    kegg_compound = get_ontology_sub_dag('br:br08001', onto_dict)
    kegg_lipids = get_ontology_sub_dag('br:br08002', onto_dict)
    kegg_peptids = get_ontology_sub_dag('br:br08005', onto_dict)
    kegg_phytochem = get_ontology_sub_dag('br:br08003', onto_dict)
    kegg_nat_toxins = get_ontology_sub_dag('br:br08009', onto_dict)
    kegg_glycosides = get_ontology_sub_dag('br:br08021', onto_dict)
    kegg_compound.update(kegg_lipids)
    kegg_compound.update(kegg_peptids)
    kegg_compound.update(kegg_phytochem)
    kegg_compound.update(kegg_nat_toxins)
    kegg_compound.update(kegg_glycosides)
    return kegg_pathway, kegg_reaction, kegg_gene_product, kegg_compound


def find_all_paths_to_root(node, parents_dict):
    if node not in parents_dict or not parents_dict[node]:
        return [[node]]
    paths = []
    for parent in parents_dict[node]:
        for path in find_all_paths_to_root(parent, parents_dict):
            paths.append([node] + path)
    return paths


INDEX = ['V-1', 'E', 'E/(V-1)', 'Density', 'Nb Leaves', 'Nb Branches', 'Max Parents',
         'Mean Parents', 'Med Parents', '%+1 Parents',
         'Nb Paths', 'Max Path len', 'Min Path len', 'Mean Path len', 'Med Path len']


def get_stats_ontology(onto_dict):
    values = pd.Series(index=INDEX)
    # Calculate total number of elements
    nb_v = len(onto_dict)
    nb_e = sum(len(v) for v in onto_dict.values())
    values['V-1'] = nb_v
    values['E'] = nb_e
    values['E/(V-1)'] = nb_e/nb_v
    values['Density'] = (2*nb_e)/(nb_v*(nb_v-1))
    # Separate leaves and branches elements
    lv_classes = set(onto_dict.keys())
    for p_classes in onto_dict.values():
        for p in p_classes:
            if p in lv_classes:
                lv_classes.remove(p)
    br_classes = set(onto_dict.keys()).difference(lv_classes)
    # Calculate number of leaves and branches elements
    values['Nb Leaves'] = len(lv_classes)
    values['Nb Branches'] = len(br_classes)
    # Calculate nb parents
    nb_parents = [len(v) for v in onto_dict.values()]
    values['Max Parents'] = np.max(nb_parents)
    # values['min_parents'] = np.min(nb_parents)
    values['Mean Parents'] = np.mean(nb_parents)
    values['Med Parents'] = np.median(nb_parents)
    # values['#>1_parents'] = len(nb_parents) - nb_parents.count(1)
    values['%+1 Parents'] = (len(nb_parents) - nb_parents.count(1))/len(nb_parents)
    # values['%=1_parents'] = nb_parents.count(1) / len(nb_parents)
    paths = list()
    for e in lv_classes:
        e_paths = find_all_paths_to_root(e, onto_dict)
        for p in e_paths:
            paths.append(len(p))
    values['Nb Paths'] = int(len(paths))
    values['Max Path len'] = int(np.max(paths))
    values['Min Path len'] = int(np.min(paths))
    values['Mean Path len'] = np.mean(paths)
    values['Med Path len'] = int(np.median(paths))
    return values, nb_parents, paths


METACYC_FILE = 'ontosunburst/Inputs/metacyc__26-0__classes.json'
KEGG_FILE = 'ontosunburst/Inputs/kegg__116-0__classes.json'
EC_FILE = 'ontosunburst/Inputs/ec__18jun25__classes.json'
CHEBI_FILE = 'ontosunburst/Inputs/chebi__239__classes.json'
CHEBI_R_FILE = 'ontosunburst/Inputs/chebi_r__239__classes.json'
GO_CC_FILE = 'ontosunburst/Inputs/go_cc__22jul25__classes.json'
GO_BP_FILE = 'ontosunburst/Inputs/go_bp__22jul25__classes.json'
GO_MF_FILE = 'ontosunburst/Inputs/go_mf__22jul25__classes.json'


def stat_all_onto():
    mc_dag = load_onto_file(METACYC_FILE)
    mc_r_dag, mc_p_dag, mc_c_dag, mc_g_dag = decompose_metacyc(mc_dag)
    kegg_dag = load_onto_file(KEGG_FILE)
    kegg_p_dag, kegg_r_dag, kegg_g_dag, kegg_c_dag = decompose_kegg(kegg_dag)
    ec_dag = load_onto_file(EC_FILE)
    chebi_dag = load_onto_file(CHEBI_FILE)
    chebi_r_dag = load_onto_file(CHEBI_R_FILE)
    go_bp_dag = load_onto_file(GO_BP_FILE)
    go_mf_dag = load_onto_file(GO_MF_FILE)
    all_dag = {'metacyc_rxn': mc_r_dag, 'metacyc_pwy': mc_p_dag, 'metacyc_cpd': mc_c_dag,
               'metacyc_gp': mc_g_dag,
               'kegg_rxn': kegg_r_dag, 'kegg_pwy': kegg_p_dag, 'kegg_cpd': kegg_c_dag,
               'kegg_gp': kegg_g_dag,
               'ec': ec_dag, 'chebi': chebi_dag, 'chebi_roles': chebi_r_dag,
               'go_bp': go_bp_dag, 'go_mf': go_mf_dag}
    # del all_dag['kegg_ko']
    # del all_dag['chebi']

    stats_df = pd.DataFrame(columns=list(all_dag.keys()), index=INDEX)
    onto_parents = {}
    onto_paths = {}
    for onto, dag in all_dag.items():
        stats, parents, paths = get_stats_ontology(dag)
        stats_df[onto] = stats
        onto_parents[onto] = parents
        onto_paths[onto] = paths

    # stats_df = stats_df.sort_values('Mean Parents', axis='columns')
    stats_df.to_csv('stats.tsv', sep='\t')
    # generate_scatter(stats_df.T, '|V|', '|E|')
    # generate_scatter(stats_df.T, '#leaves', '#branches')
    # generate_scatter(stats_df.T, '%>1_parents', 'max_parents')
    # generate_scatter(stats_df.T, 'max_path_len', 'max_parents')

    # generate_boxplot(onto_parents, 'number of parents')
    # generate_boxplot(onto_paths, 'length of paths')
    # generate_kde(onto_parents, 'number of parents')
    # generate_kde(onto_paths, 'length of paths')
    # generate_violin(onto_parents, 'number of parents')
    # generate_violin(onto_paths, 'length of paths')


def generate_boxplot(data, title):
    plt.close()
    plt.figure(figsize=(12, 6))
    plt.boxplot(data.values(), tick_labels=data.keys(), whis=1.0)
    plt.title(f'Repartition of {title}')
    y_ticks = np.arange(start=0, stop=30, step=1)
    plt.yticks(y_ticks)
    plt.xticks(rotation=45)
    plt.grid(True)
    plt.savefig(f'{title.replace(" ", "_")}_boxplot.png', dpi=600)


def generate_kde(data, title):
    plt.close()
    plt.figure(figsize=(10, 6))
    for onto, values in data.items():
        sns.kdeplot(values, label=onto, linewidth=2)
    plt.xlabel("Valeurs")
    plt.ylabel("Densité")
    plt.title(f'Repartition of {title}')
    plt.legend(title="Ontology")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'{title.replace(" ", "_")}_kde.png', dpi=600)


def generate_violin(data, title):
    plt.close()
    plt.figure(figsize=(12, 6))
    sns.violinplot(data=data)
    plt.title(f'Repartition of {title}')
    plt.xticks(rotation=45)
    plt.grid(True)
    plt.savefig(f'{title.replace(" ", "_")}_violin.png', dpi=600)


def generate_scatter(data, x, y):
    plt.close()
    plt.figure(figsize=(12, 6))
    sns.scatterplot(data=data, x=x, y=y, hue=data.index)
    plt.title(f'{x} / {y} scatter plot')
    plt.savefig(f'{x.replace(" ", "_")}_{y.replace(" ", "_")}_scatter.png', dpi=600)


stat_all_onto()