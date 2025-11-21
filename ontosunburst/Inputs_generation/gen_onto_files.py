import json
import time
import re
from kegg2bipartitegraph.reference import KEGG as KEGG_BS, extract_parent_child_nested_dict, \
    get_kegg_database_version
# from bioservices import KEGG as KEGG_BS
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
    _, k2b_dict, labels_dict = extract_kegg_metabolism_onto(labels, brite_regex)
    version = get_kegg_database_version().split('+')[0].replace('.', '-')
    sub_roots = get_sub_roots(k2b_dict)
    for sub_root in sub_roots:
        k2b_dict[sub_root] = [ROOTS[KEGG]]
    output_classes = get_output_path(KEGG, version, CLASSES_SUFFIX)
    output_labels = get_output_path(KEGG, version, LABELS_SUFFIX)
    with open(output_classes, 'w') as cf, open(output_labels, 'w') as lf:
        json.dump(fp=cf, obj=k2b_dict, indent=1)
        json.dump(fp=lf, obj=labels_dict, indent=1)


def extract_kegg_metabolism_onto(labels, brite_regex, output_file=None):
    global_hierarchy = {}
    for b_id, b_rg in brite_regex.items():
        hierarchy, labels = extract_brite_hierarchy(b_id, b_rg, labels)
        global_hierarchy[b_id] = hierarchy
    child_parent_dict = {}
    extract_parent_child_nested_dict(global_hierarchy, child_parent_dict)

    if output_file is not None:
        with open(output_file, 'w') as open_json_file:
            json.dump(global_hierarchy, open_json_file, indent=4)
        with open(output_file.replace('.json', '_labels.json'), 'w') as open_json_file:
            json.dump(labels, open_json_file, indent=4)
    return global_hierarchy, child_parent_dict, labels


def extract_brite_hierarchy(brite_id, regex, labels):
    response_text = KEGG_BIOSERVICES.get(brite_id)
    hierarchy = {}
    sep = ' -- '
    for line in response_text.split('\n'):
        # A LEVEL
        if line.startswith('A'):
            a_id = brite_id + sep + line[1:]
            a_label = line[1:]
            labels[a_id] = a_label
            hierarchy[a_id] = {}
        # B LEVEL
        if line.startswith('B  '):
            b_id = brite_id + sep + line.replace('B  ', '')
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
        # C LEVEL
        if line.startswith('C    '):
            c_id = brite_id + sep + line.replace('C    ', '')
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
        # D LEVEL
        if line.startswith('D      '):
            d_id = brite_id + sep + line.replace('D      ', '')
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
        # E LEVEL
        if line.startswith('E        '):
            e_id = brite_id + sep + line.replace('E        ', '')
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

generate_kegg_input()

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

generate_go_input(version='22jul25')