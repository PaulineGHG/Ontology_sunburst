# Ontology_sunburst
tututu
## Requirements

Requirements from `requirements.txt`

Install aucomana : https://github.com/PaulineGHG/aucomana.git
(or comment code if no need of pathway)


## Utilisation

### Inputs

#### Objects sets

Codes to extract some objects (metabolites, pathways, chebi ID, ...) in 
`obj_extraction.py` used as input for sunburst creation.

#### Ontology files

For metabolites and pathways:

- Class file : `Input/classes.json`
- Padmet ref : `Input/metacyc_26.0_prot70.padmet`

For enzymes:

- Ontology : `Input/enzymes_ontology.json`
- Names : `Input/enzymes_class_names.json`

### Run

Codes to run workflows (proportion or comparison) to create sunburst in
`class_metabolites.py`

#### Example 

```
CLASS_FILE = 'Inputs/classes.json'
P_REF = 'Inputs/metacyc_26.0_prot70.padmet'
INPUT_FILE = 'Path/To/clusters.tsv'

CLUST_CYAN = [7, 8]
CLUST_PINK = [3, 4, 5, 6]
CLUST_ALL = list(range(1, 12))

MET_CYAN = extract_metabolites_clusters(INPUT_FILE, CLUST_CYAN)
MET_PINK = extract_metabolites_clusters(INPUT_FILE, CLUST_PINK)
MET_ALL = extract_metabolites_clusters(INPUT_FILE, CLUST_ALL)

proportion_workflow(MET_CYAN, CLASS_FILE, P_REF, 'Output/prop_cyan', True)
proportion_workflow(MET_PINK, CLASS_FILE, P_REF, 'Output/prop_pink', True)
proportion_workflow(MET_ALL, CLASS_FILE, P_REF, 'Output/prop_all', True)

comparison_workflow(MET_CYAN, MET_ALL, CLASS_FILE, P_REF, 'Output/diff_cyan_all_bin', 'Binomial')
comparison_workflow(MET_PINK, MET_ALL, CLASS_FILE, P_REF, 'Output/diff_pink_all_bin', 'Binomial')
comparison_workflow(MET_CYAN, MET_ALL, CLASS_FILE, P_REF, 'Output/diff_cyan_all_hyp', 'Hypergeometric')
comparison_workflow(MET_PINK, MET_ALL, CLASS_FILE, P_REF, 'Output/diff_pink_all_hyp', 'Hypergeometric')
```