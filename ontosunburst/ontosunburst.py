from time import time
import plotly.graph_objects as go

from ontosunburst.input_preprocessing import *
from ontosunburst.dag2tree import *
from ontosunburst.onto2dag import *
# from ontosunburst.tree2sunburst import *


logging.basicConfig(level=logging.INFO)


# ==================================================================================================
#                                            WORKFLOW
# ==================================================================================================
def ontosunburst(interest: Input,
                 reference: Input | None = None,
                 ontology: OntologyName = None,
                 analysis: str = TOPOLOGY_A,
                 output: FilePath = 'sunburst',
                 scores: Dict[str, float] = None,
                 write_output: bool = True,
                 ontology_dag_input: str | OntologyDAG = None,
                 id_to_label_input: str | IdToLabel = None,
                 use_labels: bool = True,
                 test: EnrichmentTest = HYPERGEO_TEST,
                 root_cut: str = ROOT_CUT,
                 path_cut: str = PATH_UNCUT,
                 ref_base: bool = False,
                 hide_leaves: bool = False,
                 **kwargs) -> go.Figure:
    """ Main function to be called generating the sunburst figure

    Parameters
    ----------
    interest:  list[str] | set[str] | dict[str, float]
        Interest list or set of concepts IDs to classify.
        Can be associated to a weigh as a dictionary.
    reference: list[str] | set[str] | dict[str, float] | None (optional, default=None)
        Reference list or set of concepts IDs to classify.
        Can be associated to a weigh as a dictionary.
    ontology: str (optional, default=None, values in ['metacyc', 'ec', 'chebi', 'chebi_r', 'kegg',
                                                      'go_cc', 'go_bp', 'go_mf', 'go', None])
        Ontology name to use.
    analysis: str (optional, default='topology', values in ['topology', 'enrichment'])
        Analysis mode : topology or enrichment.
    output: str (optional, default='sunburst')
        Path of the output to save figure, if None, outputs will be sunburst.html and sunburst.tsv
        files
    scores: Dict[str, float] (optional, default=None)
        Dictionary associating for each ontology ID, its precalculated enrichment score. If None
        enrichment will be calculated.
    write_output: bool (optional, default=True)
        True to write the html figure and tsv class files, False to only return plotly sunburst
        figure.
    ontology_dag_input: str or Dict[str, str] (optional, default=None)
        Ontology DAG dictionary or json file. Use if tailored ontology or alternative
        (modified, updated, ...) default ontology DAG.
    id_to_label_input: str or Dict[str, str] (optional, default=None)
        Path to ID-LABELS association json file or ID-LABELS association dictionary.
        If None default files will be used. Use if tailored ontology or alternative
        (modified, updated, ...) default ontology.
    labels: bool (optional, default=True)
        True to show labels as sunburst sectors labels, False to show ID as sunburst sectors labels.
    test: str (optional, default='binomial', values in ['binomial', 'hypergeometric'])
        Type of test if analysis=enrichment, binomial or hypergeometric test.
    root_cut: str (optional, default='cut', values in ['uncut', 'cut', 'total'])
        mode for root cutting (uncut, cut or total)
    path_cut: str (optional, default='uncut', values in ['uncut', 'deeper', 'higher', 'bound'])
        mode for nested path cutting (uncut, deeper, higher or bound)
    ref_base: bool (optional, default=False)
        True to have the base classes representation of the reference set in the figure.
    show_leaves: bool (optional, default=False)
        True to show input metabolic objets at sunburst leaves
    **kwargs

    Returns
    -------
    go.Figure
        Plotly graph_objects figure of the sunburst
    """
    start_time = time()

    # ================================== INPUT PREPROCESSING =======================================
    # MANAGE INPUTS --------------------------------------------------------------------------------
    interest, reference = check_inputs_sets(interest, reference)
    all_concepts = set(interest.keys()).union(set(reference.keys()))
    # LOAD ONTOLOGY DAG DICTIONARY -----------------------------------------------------------------
    ontology_dag = get_ontology_dag_dict(ontology, ontology_dag_input)
    detect_cycles(ontology_dag)
    unmapped = check_input_ids_to_ontology_mapping(ontology_dag, all_concepts)
    interest, reference, all_concepts = manage_unmapped(interest, reference, unmapped)
    # GET ROOT -------------------------------------------------------------------------------------
    root, ontology_dag = get_ontology_root(ontology_dag)
    # LOAD ID TO LABELS DICTIONARY -----------------------------------------------------------------
    id_to_label = get_id_to_label_dict(ontology, id_to_label_input)

    # ===================================== ONTO TO DAG ============================================
    # GENERATE SUB-DAG FROM INPUT ------------------------------------------------------------------
    sub_dag = SubDAG(interest, reference, all_concepts, ontology_dag, root, id_to_label, test)

    end_time = time()
    logging.info(f'Execution time : {end_time - start_time} seconds')
    # return fig


# ==================================================================================================
#                                             FUNCTIONS
# ==================================================================================================
def write_concepts_classes(ontology: str, all_classes: Dict[str, Set[str]], output: str,
                           id_to_label: Dict[str, str]):
    """ Writes, for each input class, all its ancestors in a .tsv file.

    Parameters
    ----------
    ontology
    all_classes
    output
    id_to_label
    """
    if ontology is None:
        ontology = ''
    links_dict = {METACYC: 'https://metacyc.org/compound?orgid=META&id=',
                  CHEBI: 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:',
                  CHEBI_R: 'https://www.ebi.ac.uk/chebi/searchId.do?chebiId=CHEBI:',
                  EC: 'https://enzyme.expasy.org/EC/',
                  KEGG: 'https://www.genome.jp/entry/',
                  GO_MF: 'https://amigo.geneontology.org/amigo/term/',
                  GO_CC: 'https://amigo.geneontology.org/amigo/term/',
                  GO_BP: 'https://amigo.geneontology.org/amigo/term/',
                  GO: 'https://amigo.geneontology.org/amigo/term/',
                  '': ''}
    with open(f'{output}.tsv', 'w') as f:
        f.write('\t'.join(['ID', 'Label', 'Classes ID', 'Classes Label', 'Link']) + '\n')
        for met_id, classes_id, in all_classes.items():
            link = links_dict[ontology] + met_id
            met_lab = get_name(met_id, id_to_label)
            classes_lab = [get_name(cl, id_to_label) for cl in classes_id]
            f.write('\t'.join([met_id, met_lab, ', '.join(classes_id), ', '.join(classes_lab),
                               link]) + '\n')