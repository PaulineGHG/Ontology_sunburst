import json
import time
import re
import kegg2bipartitegraph.reference as keggr
from bioservices import KEGG as KEGG_BS
from SPARQLWrapper import SPARQLWrapper, JSON
from padmet.classes.padmetRef import PadmetRef
from ontosunburst.ontosunburst import METACYC, EC, CHEBI, CHEBI_R, GO_MF, GO_BP, GO_CC, KEGG, \
    ROOTS, CLASSES_SUFFIX, LABELS_SUFFIX


KEGG_BIOSERVICES = KEGG_BS()
URL = {CHEBI: 'http://localhost:3030/chebi/',
       GO_MF: 'http://localhost:3030/go/'}

GO_ROOTS = [GO_CC, GO_BP, GO_MF]


def get_output_path(prefix, version, suffix):
    return prefix + '__' + version + '__' + suffix


def get_sub_roots(dict_onto):
    sub_roots = set()
    for v in dict_onto.values():
        for c in v:
            if c not in dict_onto:
                sub_roots.add(c)
    return sub_roots


# METACYC
# ==================================================================================================
def generate_metacyc_input(input_padmet, version=''):
    output = get_output_path(METACYC, version, CLASSES_SUFFIX)
    pref = PadmetRef(input_padmet)
    rels = pref.getAllRelation()
    mc_classes = dict()
    for r in rels:
        if r.type == 'is_a_class':
            if r.id_in not in mc_classes:
                mc_classes[r.id_in] = set()
            mc_classes[r.id_in].add(r.id_out)
    for c, p in mc_classes.items():
        mc_classes[c] = list(p)
    with open(output, 'w') as o:
        json.dump(fp=o, obj=mc_classes, indent=1)


# EC
# ==================================================================================================
def generate_ec_input(enzclass_txt, enzyme_dat, version=''):
    output_name = get_output_path(EC, version, LABELS_SUFFIX)
    output_class = get_output_path(EC, version, CLASSES_SUFFIX)
    names = dict()
    classes = dict()
    with open(enzclass_txt, 'r') as f:
        for l in f:
            ec_id = l[:9].replace(' ', '')
            if ec_id.count('.') == 3:
                ec_name = l[10:].strip()
                names[ec_id] = ec_name
                ec_parent = get_ec_parent(ec_id)
                classes[ec_id] = [ec_parent]
    with open(enzyme_dat, 'r') as f:
        for l in f:
            if l.startswith('ID'):
                ec_id = l[3:].strip()
            if l.startswith('DE'):
                ec_name = l[3:].strip()
                names[ec_id] = ec_name
                ec_parent = get_ec_parent(ec_id)
                classes[ec_id] = [ec_parent]
    with open(output_name, 'w') as on, open(output_class, 'w') as oc:
        json.dump(fp=on, obj=names, indent=1)
        json.dump(fp=oc, obj=classes, indent=1)


def get_ec_parent(ec: str) -> str:
    ec_lvl = ec.count('-')
    if ec_lvl == 3:
        return ROOTS[EC]
    else:
        ec_lst = ec.split('.')
        ec_lst[-(ec_lvl+1)] = '-'
        return '.'.join(ec_lst)


# KEGG
# ==================================================================================================
def generate_kegg_input():
    keggr.create_reference_base()
    print('base created')
    _, k2b_dict = keggr.get_kegg_hierarchy()
    version = keggr.get_kegg_database_version().split('+')[0].replace('.', '-')
    sub_roots = get_sub_roots(k2b_dict)
    for sub_root in sub_roots:
        k2b_dict[sub_root] = [ROOTS[KEGG]]
    output = get_output_path(KEGG, version, CLASSES_SUFFIX)
    with open(output, 'w') as f:
        json.dump(fp=f, obj=k2b_dict, indent=1)

    kegg_hierarchy = {'br:ko00002': 'Pathways',
                      'br:br08901': 'Pathways',
                      'br:ko00001': 'Gene products',
                      'br:br08005': 'Gene products',
                      'br:br08001': 'Compounds',
                      'br:br08002': 'Compounds',
                      'br:br08003': 'Compounds',
                      'br:br08009': 'Compounds',
                      'br:br08021': 'Compounds',
                      'br:br08201': 'Reactions',
                      'Pathways': 'KEGG',
                      'Gene products': 'KEGG',
                      'Compounds': 'KEGG',
                      'Reactions': 'KEGG'}


def extract_kegg_metabolism_onto(output_file):
    labels = {'br:ko00002': 'KEGG modules',
              'br:br08901': 'KEGG pathway maps',
              'br:ko00001': 'KEGG Orthology (KO)',
              'br:br08005': 'Bioactive peptides',
              'br:br08001': 'Compounds with biological roles',
              'br:br08002': 'Lipids',
              'br:br08003': 'Phytochemical compounds',
              'br:br08009': 'Natural toxins',
              'br:br08021': 'Glycosides',
              'br:br08201': 'Enzymatic reactions'}

    brite_regex = {'br:ko00002': r'M\d{5}',
                   'br:br08901': r'\d{5}',
                   'br:ko00001': r'K\d{5}',
                   'br:br08005': r'C\d{5}',
                   'br:br08001': r'C\d{5}',
                   'br:br08002': r'[GC]\d{5}',
                   'br:br08003': r'C\d{5}',
                   'br:br08009': r'C\d{5}',
                   'br:br08021': r'C\d{5}',
                   'br:br08201': r'R\d{5}'}

    global_hierarchy = {}
    for b_id, b_rg in brite_regex.items():
        hierarchy, labels = extract_brite_hierarchy(b_id, b_rg, labels)
        global_hierarchy[b_id] = hierarchy

    if output_file is not None:
        with open(output_file, 'w') as open_json_file:
            json.dump(global_hierarchy, open_json_file, indent=4)
        with open(output_file.replace('.json', '_labels.json'), 'w') as open_json_file:
            json.dump(labels, open_json_file, indent=4)


def extract_brite_hierarchy(brite_id, regex, labels):
    response_text = KEGG_BIOSERVICES.get(brite_id)
    hierarchy = {}
    for line in response_text.split('\n'):
        if line.startswith('A'):
            a_id = brite_id + '_' + line[1:]
            a_label = line[1:]
            labels[a_id] = a_label
            hierarchy[a_id] = {}
        if line.startswith('B  '):
            b_id = brite_id + '_' + line.replace('B  ', '')
            b_label = line.replace('B  ', '')
            if b_id == a_id:
                b_id = b_id + ' B'
            e_b_id = b_label.split('  ')[0]
            if re.match(regex, e_b_id):
                e_b_label = b_label.split('  ')[1]
                labels[e_b_id] = e_b_label
                if hierarchy[a_id] == {}:
                    hierarchy[a_id] = []
                hierarchy[a_id].append(e_b_id)
            else:
                labels[b_id] = b_label
                hierarchy[a_id][b_id] = {}
        if line.startswith('C    '):
            c_id = brite_id + '_' + line.replace('C    ', '')
            c_label = line.replace('C    ', '')
            if c_id in [a_id, b_id]:
                c_id = c_id + ' C'
            e_c_id = c_label.split('  ')[0]
            if re.match(regex, e_c_id):
                e_c_label = c_label.split('  ')[1]
                labels[e_c_id] = e_c_label
                if hierarchy[a_id][b_id] == {}:
                    hierarchy[a_id][b_id] = []
                hierarchy[a_id][b_id].append(e_c_id)
            else:
                labels[c_id] = c_label
                hierarchy[a_id][b_id][c_id] = {}
        if line.startswith('D      '):
            d_id = brite_id + '_' + line.replace('D      ', '')
            d_label = line.replace('D      ', '')
            if d_id in [a_id, b_id, c_id]:
                d_id = d_id + ' D'
            e_d_id = d_label.split('  ')[0]
            if re.match(regex, e_d_id):
                e_d_label = d_label.split('  ')[1]
                labels[e_d_id] = e_d_label
                if hierarchy[a_id][b_id][c_id] == {}:
                    hierarchy[a_id][b_id][c_id] = []
                hierarchy[a_id][b_id][c_id].append(e_d_id)
            else:
                labels[d_id] = d_label
                hierarchy[a_id][b_id][c_id][d_id] = {}
        if line.startswith('E        '):
            e_id = brite_id + '_' + line.replace('E        ', '')
            e_label = line.replace('E        ', '')
            if e_id in [a_id, b_id, c_id, d_id]:
                e_id = e_id + ' E'
            e_e_id = e_label.split('  ')[0]
            if re.match(regex, e_e_id):
                e_e_label = e_label.split('  ')[1]
                labels[e_e_id] = e_e_label
                if hierarchy[a_id][b_id][c_id][d_id] == {}:
                    hierarchy[a_id][b_id][c_id][d_id] = []
                hierarchy[a_id][b_id][c_id][d_id].append(e_e_id)
            else:
                '|!!!| max depth reached'
    time.sleep(3)
    return hierarchy, labels


def get_kegg_hierarchy(hierarchy_file=None):
    """Using bioservices.KEGG to create hierarchy of KEGG IDs.

    Args:
        hierarchy_file (str): output file which will contains hierarchy of KEGG IDs
    Returns:
        global_hierarchy (dict): nested dictionary having hierarchy of KEGG IDs from brite hierarchy (higher elements
        in dictionary correspond to higher element in hierarchy)
        parent_child_dict (dict): dictionary with child KEGG ID as key and a list of their direct parent KEGG elements
        as value
    """
    # PWY AND MODULES
    # ----------------------------------------------------------
    # Create hierarchy for module ID.
    response_text = KEGG_BIOSERVICES.get('br:ko00002')

    module_hierarchy = {}
    for line in response_text.split('\n'):
        if line.startswith('A'):
            a_category = line[1:]
            module_hierarchy[a_category] = {}
        if line.startswith('B  '):
            b_category = line.replace('B  ', '')
            # If statement to avoid having the same module name in both parent and child.
            if b_category == a_category:
                b_category = b_category + ' B'
            module_hierarchy[a_category][b_category] = {}
        if line.startswith('C    '):
            c_category = line.replace('C    ', '')
            if c_category in [a_category, b_category]:
                c_category = c_category + ' C'
            module_hierarchy[a_category][b_category][c_category] = []
        if line.startswith('D      '):
            d_category = line.replace('D      ', '')
            module_id = d_category.split('  ')[0]
            module_hierarchy[a_category][b_category][c_category].append(module_id)
    time.sleep(3)

    # Create hierarchy for pathway.
    response_text = KEGG_BIOSERVICES.get('br:br08901')
    pathway_hierarchy = {}
    for line in response_text.split('\n'):
        if line.startswith('A'):
            a_category = line[1:]
            pathway_hierarchy[a_category] = {}
        if line.startswith('B  '):
            b_category = line.replace('B  ', '')
            if b_category == a_category:
                b_category = b_category + ' B'
            pathway_hierarchy[a_category][b_category] = []
        if line.startswith('C    '):
            pathway_id = 'map' + line.replace('C    ', '').split('  ')[0]
            pathway_hierarchy[a_category][b_category].append(pathway_id)
    time.sleep(3)

    # GENE PRODUCTS
    # ----------------------------------------------------------
    # Create hierarchy for KO hierarchy.
    response_text = KEGG_BIOSERVICES.get('br:ko00001')
    ko_hierarchy = {}
    for line in response_text.split('\n'):
        if line.startswith('A'):
            a_category = line[1:]
            ko_hierarchy[a_category] = {}
        if line.startswith('B  '):
            b_category = line.replace('B  ', '')
            if b_category == a_category:
                b_category = b_category + ' B'
            ko_hierarchy[a_category][b_category] = {}
        if line.startswith('C    '):
            c_category = line.replace('C    ', '')
            if c_category in [a_category, b_category]:
                c_category = c_category + ' C'
            ko_hierarchy[a_category][b_category][c_category] = []
        if line.startswith('D      '):
            d_category = line.replace('D      ', '')
            ko_id = d_category.split('  ')[0]
            ko_hierarchy[a_category][b_category][c_category].append(ko_id)
    time.sleep(3)

    # Create hierarchy for Bioactive peptides.
    response_text = KEGG_BIOSERVICES.get('br:br08005')
    metabolite_regex = r'C\d{5}'
    peptide_hierarchy = {}
    for line in response_text.split('\n'):
        if line.startswith('A'):
            a_category = line[1:]
            peptide_hierarchy[a_category] = {}
        if line.startswith('B  '):
            b_category = line.replace('B  ', '')
            if b_category == a_category:
                b_category = b_category + ' B'
            peptide_hierarchy[a_category][b_category] = {}
        if line.startswith('C    '):
            c_category = line.replace('C    ', '')
            if c_category in [a_category, b_category]:
                c_category = c_category + ' C'
            peptide_c_category = c_category.split('  ')[0]
            if re.match(metabolite_regex, peptide_c_category):
                if peptide_hierarchy[a_category][b_category] == {}:
                    peptide_hierarchy[a_category][b_category] = []
                peptide_hierarchy[a_category][b_category].append(peptide_c_category)
            else:
                peptide_hierarchy[a_category][b_category][c_category] = []
        if line.startswith('D      '):
            d_category = line.replace('D      ', '')
            metabolite_id = d_category.split('  ')[0]
            peptide_hierarchy[a_category][b_category][c_category].append(metabolite_id)
    time.sleep(3)

    # COMPOUNDS
    # ----------------------------------------------------------
    # Create hierarchy for metabolite.
    response_text = KEGG_BIOSERVICES.get('br:br08001')
    metabolite_regex = r'C\d{5}'
    metabolite_hierarchy = {}
    for line in response_text.split('\n'):
        if line.startswith('A'):
            a_category = line[1:]
            metabolite_hierarchy[a_category] = {}
        if line.startswith('B  '):
            b_category = line.replace('B  ', '')
            if b_category == a_category:
                b_category = b_category + ' B'
            metabolite_hierarchy[a_category][b_category] = {}
        if line.startswith('C    '):
            c_category = line.replace('C    ', '')
            if c_category in [a_category, b_category]:
                c_category = c_category + ' C'
            metabolite_c_category = c_category.split('  ')[0]
            if re.match(metabolite_regex, metabolite_c_category):
                if metabolite_hierarchy[a_category][b_category] == {}:
                    metabolite_hierarchy[a_category][b_category] = []
                metabolite_hierarchy[a_category][b_category].append(metabolite_c_category)
            else:
                metabolite_hierarchy[a_category][b_category][c_category] = []
            metabolite_hierarchy[a_category][b_category][c_category] = []
        if line.startswith('D      '):
            d_category = line.replace('D      ', '')
            metabolite_id = d_category.split('  ')[0]
            metabolite_hierarchy[a_category][b_category][c_category].append(metabolite_id)
    time.sleep(3)

    # Create hierarchy for lipids.
    response_text = KEGG_BIOSERVICES.get('br:br08002')
    lipid_regex = r'[GC]\d{5}'
    lipid_hierarchy = {}
    for line in response_text.split('\n'):
        if line.startswith('A'):
            a_category = line[1:]
            lipid_hierarchy[a_category] = {}
        if line.startswith('B  '):
            b_category = line.replace('B  ', '')
            if b_category == a_category:
                b_category = b_category + ' B'
            lipid_hierarchy[a_category][b_category] = {}
        if line.startswith('C    '):
            c_category = line.replace('C    ', '')
            if c_category in [a_category, b_category]:
                c_category = c_category + ' C'
            lipid_c_category = c_category.split('  ')[0]
            if re.match(lipid_regex, lipid_c_category):
                if lipid_hierarchy[a_category][b_category] == {}:
                    lipid_hierarchy[a_category][b_category] = []
                lipid_hierarchy[a_category][b_category].append(lipid_c_category)
            else:
                lipid_hierarchy[a_category][b_category][c_category] = {}
        if line.startswith('D      '):
            d_category = line.replace('D      ', '')
            if d_category in [a_category, b_category, c_category]:
                d_category = d_category + ' D'
            lipid_id = d_category.split('  ')[0]
            if re.match(lipid_regex, lipid_id):
                if lipid_hierarchy[a_category][b_category][c_category] == {}:
                    lipid_hierarchy[a_category][b_category][c_category] = []
                lipid_hierarchy[a_category][b_category][c_category].append(lipid_id)
            else:
                lipid_hierarchy[a_category][b_category][c_category][d_category] = []
        if line.startswith('E        '):
            e_category = line.replace('E        ', '')
            lipid_id = e_category.split('  ')[0]
            if re.match(lipid_regex, lipid_id):
                lipid_hierarchy[a_category][b_category][c_category][d_category].append(lipid_id)
    time.sleep(3)

    # Create hierarchy for Phytochemical compounds.
    response_text = KEGG_BIOSERVICES.get('br:br08003')
    phytocpd_hierarchy = {}
    for line in response_text.split('\n'):
        if line.startswith('A'):
            a_category = line[1:]
            phytocpd_hierarchy[a_category] = {}
        if line.startswith('B  '):
            b_category = line.replace('B  ', '')
            if b_category == a_category:
                b_category = b_category + ' B'
            phytocpd_hierarchy[a_category][b_category] = {}
        if line.startswith('C    '):
            c_category = line.replace('C    ', '')
            if c_category in [a_category, b_category]:
                c_category = c_category + ' C'
            metabolite_c_category = c_category.split('  ')[0]
            if re.match(metabolite_regex, metabolite_c_category):
                if phytocpd_hierarchy[a_category][b_category] == {}:
                    phytocpd_hierarchy[a_category][b_category] = []
                phytocpd_hierarchy[a_category][b_category].append(metabolite_c_category)
            else:
                phytocpd_hierarchy[a_category][b_category][c_category] = []
        if line.startswith('D      '):
            d_category = line.replace('D      ', '')
            metabolite_id = d_category.split('  ')[0]
            phytocpd_hierarchy[a_category][b_category][c_category].append(metabolite_id)
    time.sleep(3)

    # Create hierarchy for Natural toxins.
    response_text = KEGG_BIOSERVICES.get('br:br08009')
    natox_hierarchy = {}
    for line in response_text.split('\n'):
        if line.startswith('A'):
            a_category = line[1:]
            natox_hierarchy[a_category] = {}
        if line.startswith('B  '):
            b_category = line.replace('B  ', '')
            if b_category == a_category:
                b_category = b_category + ' B'
            natox_hierarchy[a_category][b_category] = {}
        if line.startswith('C    '):
            c_category = line.replace('C    ', '')
            if c_category in [a_category, b_category]:
                c_category = c_category + ' C'
            metabolite_c_category = c_category.split('  ')[0]
            if re.match(metabolite_regex, metabolite_c_category):
                if natox_hierarchy[a_category][b_category] == {}:
                    natox_hierarchy[a_category][b_category] = []
                natox_hierarchy[a_category][b_category].append(metabolite_c_category)
            else:
                natox_hierarchy[a_category][b_category][c_category] = []
        if line.startswith('D      '):
            d_category = line.replace('D      ', '')
            metabolite_id = d_category.split('  ')[0]
            natox_hierarchy[a_category][b_category][c_category].append(metabolite_id)
    time.sleep(3)

    # Create hierarchy for Glycosides.
    response_text = KEGG_BIOSERVICES.get('br:br08021')
    glycoside_hierarchy = {}
    for line in response_text.split('\n'):
        if line.startswith('A'):
            a_category = line[1:]
            glycoside_hierarchy[a_category] = {}
        if line.startswith('B  '):
            b_category = line.replace('B  ', '')
            if b_category == a_category:
                b_category = b_category + ' B'
            glycoside_hierarchy[a_category][b_category] = {}
        if line.startswith('C    '):
            c_category = line.replace('C    ', '')
            if c_category in [a_category, b_category]:
                c_category = c_category + ' C'
            metabolite_c_category = c_category.split('  ')[0]
            if re.match(metabolite_regex, metabolite_c_category):
                if glycoside_hierarchy[a_category][b_category] == {}:
                    glycoside_hierarchy[a_category][b_category] = []
                glycoside_hierarchy[a_category][b_category].append(metabolite_c_category)
            else:
                glycoside_hierarchy[a_category][b_category][c_category] = []
        if line.startswith('D      '):
            d_category = line.replace('D      ', '')
            metabolite_id = d_category.split('  ')[0]
            glycoside_hierarchy[a_category][b_category][c_category].append(metabolite_id)
    time.sleep(3)

    # RXN
    # ----------------------------------------------------------
    # Create hierarchy for reactions.
    response_text = KEGG_BIOSERVICES.get('br:br08201')
    rxn_regex = r'R\d{5}'
    rxn_hierarchy = {}
    for line in response_text.split('\n'):
        if line.startswith('A'):
            a_category = line[1:]
            rxn_hierarchy[a_category] = {}
        if line.startswith('B  '):
            b_category = line.replace('B  ', '')
            if b_category == a_category:
                b_category = b_category + ' B'
            rxn_hierarchy[a_category][b_category] = {}
        if line.startswith('C    '):
            c_category = line.replace('C    ', '')
            if c_category in [a_category, b_category]:
                c_category = c_category + ' C'
            lipid_c_category = c_category.split('  ')[0]
            if re.match(rxn_regex, lipid_c_category):
                if rxn_hierarchy[a_category][b_category] == {}:
                    rxn_hierarchy[a_category][b_category] = []
                rxn_hierarchy[a_category][b_category].append(lipid_c_category)
            else:
                rxn_hierarchy[a_category][b_category][c_category] = {}
        if line.startswith('D      '):
            d_category = line.replace('D      ', '')
            if d_category in [a_category, b_category, c_category]:
                d_category = d_category + ' D'
            lipid_id = d_category.split('  ')[0]
            if re.match(rxn_regex, lipid_id):
                if rxn_hierarchy[a_category][b_category][c_category] == {}:
                    rxn_hierarchy[a_category][b_category][c_category] = []
                rxn_hierarchy[a_category][b_category][c_category].append(lipid_id)
            else:
                rxn_hierarchy[a_category][b_category][c_category][d_category] = []
        if line.startswith('E        '):
            e_category = line.replace('E        ', '')
            lipid_id = e_category.split('  ')[0]
            if re.match(rxn_regex, lipid_id):
                rxn_hierarchy[a_category][b_category][c_category][d_category].append(lipid_id)
    time.sleep(3)

    global_hierarchy = {}
    global_hierarchy['module'] = module_hierarchy
    global_hierarchy['pathway'] = pathway_hierarchy
    global_hierarchy['ko'] = ko_hierarchy
    global_hierarchy['bioactive_peptides'] = peptide_hierarchy
    global_hierarchy['metabolite'] = metabolite_hierarchy
    global_hierarchy['lipid'] = lipid_hierarchy
    global_hierarchy['phytochemical_cpd'] = phytocpd_hierarchy
    global_hierarchy['natural_toxin'] = natox_hierarchy
    global_hierarchy['glycoside'] = glycoside_hierarchy
    global_hierarchy['enz_rxn'] = rxn_hierarchy

    if hierarchy_file is not None:
        with open(hierarchy_file, 'w') as open_json_file:
            json.dump(global_hierarchy, open_json_file, indent=4)
    return global_hierarchy

extract_kegg_metabolism_onto('truc.json')




# CHEBI CLASS
# ==================================================================================================

def chebi_onto_query(endpoint_url: str) -> dict:
    """ Returns the query results to get the chebi classes ontology tree from the role root
    (Root_id = 24431)

    Parameters
    ----------
    endpoint_url: str
        Endpoint URL of Jena Fuseki server

    Returns
    -------
    dict
        Dictionary of query results
    """
    query = f"""
            PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
            PREFIX chebidb: <http://purl.obolibrary.org/obo/CHEBI_>
            PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>

            SELECT DISTINCT ?childLabel ?parentLabel ?childId ?parentId
            WHERE {{
                VALUES ?root{{chebidb:24431}}                                        
                ?root rdfs:label ?rootLabel .

                ?child rdfs:subClassOf* ?root .
                ?child oboInOwl:id ?childId .  
                ?child rdfs:label ?childLabel .

                ?child rdfs:subClassOf ?parent .
                ?parent rdfs:label ?parentLabel .
                ?parent oboInOwl:id ?parentId . 
            }}
            """
    sparql = SPARQLWrapper(endpoint_url)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    return results


def generate_chebi_input(version='', url_endpoint=URL[CHEBI]):
    d_ontology = dict()
    d_labels = dict()
    results = chebi_onto_query(url_endpoint)
    for result in results['results']['bindings']:
        child_label = result['childLabel']['value']
        parent_label = result['parentLabel']['value']
        child_id = result['childId']['value']
        parent_id = result['parentId']['value']
        if child_id not in d_ontology:
            d_ontology[child_id] = []
        d_ontology[child_id].append(parent_id)
        d_labels[parent_id] = parent_label
        d_labels[child_id] = child_label
    sub_roots = get_sub_roots(d_ontology)
    for sub_root in sub_roots:
        d_ontology[sub_root] = [ROOTS[CHEBI]]
    output_classes = get_output_path(CHEBI, version, CLASSES_SUFFIX)
    output_labels = get_output_path(CHEBI, version, LABELS_SUFFIX)
    with open(output_classes, 'w') as oc, open(output_labels, 'w') as ol:
        json.dump(d_ontology, oc, indent=1)
        json.dump(d_labels, ol, indent=1)


# CHEBI ROLES
# ==================================================================================================
def chebi_role_onto_query(endpoint_url: str) -> dict:
    """ Returns the query results to get the chebi roles ontology tree from the role root
    (Root_id = 50906)

    Parameters
    ----------
    endpoint_url: str
        Endpoint URL of Jena Fuseki server

    Returns
    -------
    dict
        Dictionary of query results
    """
    query = f"""
            PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
            PREFIX chebidb: <http://purl.obolibrary.org/obo/CHEBI_>
            PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>

            SELECT DISTINCT ?childRoleLabel ?parentRoleLabel ?childRoleId ?parentRoleId
            WHERE {{
                VALUES ?rootRole{{chebidb:50906}}                                        
                ?rootRole rdfs:label ?rootRoleLabel .
                
                ?childRole rdfs:subClassOf* ?rootRole .
                ?childRole oboInOwl:id ?childRoleId .  
                ?childRole rdfs:label ?childRoleLabel .
                
                ?childRole rdfs:subClassOf ?parentRole .
                ?parentRole rdfs:label ?parentRoleLabel .
                ?parentRole oboInOwl:id ?parentRoleId . 
            }}
            """
    sparql = SPARQLWrapper(endpoint_url)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    return results


def chebi_chem_roles_query(endpoint_url: str) -> dict:
    """ Returns the query results to get the chebi roles associated to each chemical having a role.

    Parameters
    ----------
    endpoint_url: str
        Endpoint URL of Jena Fuseki server

    Returns
    -------
    dict
        Dictionary of query results
    """
    query = f"""
            PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
            PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
            PREFIX owl: <http://www.w3.org/2002/07/owl#>
            PREFIX obo: <http://purl.obolibrary.org/obo/>
            PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>
         
            SELECT DISTINCT ?chemLabel ?roleLabel ?chemId ?roleId
            WHERE {{                                      
                ?chem rdfs:subClassOf ?restriction .
                ?restriction rdf:type owl:Restriction .
                ?restriction owl:onProperty obo:RO_0000087 .
                ?restriction owl:someValuesFrom/(rdfs:subClassOf) ?role .
                ?role rdfs:label ?roleLabel .
                ?role oboInOwl:id ?roleId .
                ?chem rdfs:label ?chemLabel .
                ?chem oboInOwl:id ?chemId .
            }}
            """
    sparql = SPARQLWrapper(endpoint_url)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    return results


def generate_chebi_roles_input(version='', url_endpoint=URL[CHEBI]):
    d_roles_ontology = dict()
    d_labels = dict()
    results = chebi_role_onto_query(url_endpoint)
    for result in results['results']['bindings']:
        child_role_label = result['childRoleLabel']['value']
        parent_role_label = result['parentRoleLabel']['value']
        child_role_id = result['childRoleId']['value']
        parent_role_id = result['parentRoleId']['value']
        if child_role_id not in d_roles_ontology:
            d_roles_ontology[child_role_id] = []
        d_roles_ontology[child_role_id].append(parent_role_id)
        d_labels[parent_role_id] = parent_role_label
        d_labels[child_role_id] = child_role_label

    results = chebi_chem_roles_query(url_endpoint)
    for result in results['results']['bindings']:
        chem_label = result['chemLabel']['value']
        role_label = result['roleLabel']['value']
        chem_id = result['chemId']['value']
        role_id = result['roleId']['value']
        if role_id not in d_labels:
            print(f'Role {role_label} not in the ontology.')
        if chem_id not in d_roles_ontology:
            d_roles_ontology[chem_id] = []
        d_roles_ontology[chem_id].append(role_id)
        d_labels[chem_id] = chem_label

    output_classes = get_output_path(CHEBI_R, version, CLASSES_SUFFIX)
    output_labels = get_output_path(CHEBI_R, version, LABELS_SUFFIX)
    with open(output_classes, 'w') as oc, open(output_labels, 'w') as ol:
        json.dump(d_roles_ontology, oc, indent=1)
        json.dump(d_labels, ol, indent=1)


# GO
# ==================================================================================================
def go_onto_query(root: str, endpoint_url: str) -> dict:
    """ Returns the query results to get the chebi roles associated to each chemical having a role.

    Parameters
    ----------
    root: str
    endpoint_url: str
        Endpoint URL of Jena Fuseki server

    Returns
    -------
    dict
        Dictionary of query results
    """
    query = f"""
            PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
            PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
            PREFIX owl: <http://www.w3.org/2002/07/owl#>
            PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
            PREFIX dc: <http://purl.org/dc/elements/1.1/>
            PREFIX dcterms: <http://purl.org/dc/terms/>
            PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>
            PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
            PREFIX uniprot: <http://purl.uniprot.org/uniprot/>
            PREFIX up:<http://purl.uniprot.org/core/>
            PREFIX go: <http://purl.obolibrary.org/obo/GO_>
            PREFIX goavoc: <http://bio2rdf.org/goa_vocabulary:>
    
            SELECT ?childLabel ?parentLabel ?childId ?parentId
            WHERE {{
                VALUES ?root{{{root.lower()}}}
                ?root rdfs:label ?rootLabel .
                
                ?child rdfs:subClassOf* ?root .
                ?child oboInOwl:id ?childId .  
                ?child rdfs:label ?childLabel .

                ?child rdfs:subClassOf ?parent .
                ?parent rdfs:label ?parentLabel .
                ?parent oboInOwl:id ?parentId .
            }}
            """
    sparql = SPARQLWrapper(endpoint_url)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    return results


def generate_go_input(version='', url_endpoint=URL[GO_MF]):
    for root_name in GO_ROOTS:
        root = ROOTS[root_name]
        d_ontology = dict()
        d_labels = dict()
        results = go_onto_query(root, url_endpoint)
        for result in results['results']['bindings']:
            child_label = result['childLabel']['value']
            parent_label = result['parentLabel']['value']
            child_id = result['childId']['value']
            parent_id = result['parentId']['value']
            if child_id not in d_ontology:
                d_ontology[child_id] = []
            d_ontology[child_id].append(parent_id)
            d_labels[parent_id] = parent_label
            d_labels[child_id] = child_label

        output_classes = get_output_path(root_name, version, CLASSES_SUFFIX)
        output_labels = get_output_path(root_name, version, LABELS_SUFFIX)
        with open(output_classes, 'w') as oc, open(output_labels, 'w') as ol:
            json.dump(d_ontology, oc, indent=1)
            json.dump(d_labels, ol, indent=1)
