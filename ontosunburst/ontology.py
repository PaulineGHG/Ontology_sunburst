from typing import List, Set, Dict, Tuple
from SPARQLWrapper import SPARQLWrapper, JSON


# CONSTANTS ========================================================================================


METACYC = 'metacyc'
EC = 'ec'
CHEBI = 'chebi'
GO = 'go'
KEGG = 'kegg'

ROOTS = {METACYC: 'FRAMES',
         CHEBI: 'role',
         EC: 'Enzyme',
         GO: 'GO',
         KEGG: 'kegg'}

GO_ROOTS = ['cellular_component', 'biological_process', 'molecular_function']


# For MetaCyc and Kegg Ontology
# --------------------------------------------------------------------------------------------------
def extract_met_classes(metabolic_objects: List[str],
                        d_classes_ontology: Dict[str, List[str]]) -> Dict[str, List[str]]:
    """ Extract +1 parent classes for each metabolic object.

    Parameters
    ----------
    metabolic_objects: List[str]
        List of metabolic objects.
    d_classes_ontology: Dict[str, List[str]]
        Ontology classes dictionary : Dict[object, List[parents]]

    Returns
    -------
    Dict[str, List[str]], (Dict[object, List[parents]])
        Dictionary associating for each metabolic object, the list of +1 parent classes it belongs
        to.
    """
    d_obj_classes = dict()
    print(f'{len(metabolic_objects)} metabolic objects to classify')
    for obj in metabolic_objects:
        try:
            d_obj_classes[obj] = d_classes_ontology[obj]
            classified = True
        except KeyError:
            classified = False
        if not classified:
            print(f'{obj} not classified.')
    print(f'{len(d_obj_classes)}/{len(metabolic_objects)} metabolic objects classified')
    return d_obj_classes


# For EC
# --------------------------------------------------------------------------------------------------
def extract_ec_classes(ec_list: List[str], d_classes_ontology: Dict[str, List[str]]) \
        -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """ Extract +1 parent classes for each EC number.

    Parameters
    ----------
    ec_list: List[str]
        List of EC numbers
    d_classes_ontology: Dict[str, List[str]]
        EC Ontology classes dictionary : Dict[object, List[parents]]

    Returns
    -------
    Dict[str, List[str]]
        Dictionary associating for each metabolic object, the list of +1 parent classes it belongs
        to.
    Dict[str, List[str]]
        EC Ontology classes dictionary : Dict[object, List[parents]]
    """
    print(f'{len(ec_list)} EC numbers to classify')
    ec_classes = dict()
    for ec in ec_list:
        dec_ec = ec.split('.')
        while len(dec_ec) < 4:
            dec_ec.append('-')
        if dec_ec.count('-') == 3 and ec in d_classes_ontology:
            parent = ROOTS[EC]
        else:
            for i in range(3):
                if dec_ec.count('-') == i:
                    dec_ec[3-i] = '-'
                    break
            parent = '.'.join(dec_ec)
        if parent in d_classes_ontology or parent == ROOTS[EC]:
            d_classes_ontology[ec] = [parent]
            ec_classes[ec] = [parent]
        else:
            print(f'{ec} not classified')
    print(f'{len(ec_classes)}/{len(ec_list)} EC numbers classified')
    return ec_classes, d_classes_ontology


# For ChEBI Ontology
# --------------------------------------------------------------------------------------------------
def extract_chebi_roles(chebi_ids: List[str], endpoint_url: str) \
        -> Tuple[Dict[str, Set[str]], Dict[str, List[str]]]:
    """

    Parameters
    ----------
    chebi_ids: List[str]
        List of ChEBI IDs to extract roles associated
    endpoint_url: str
        URL endpoint string

    Returns
    -------
    Tuple[Dict[str, Set[str]], Dict[str, List[str]]]
        Dictionary of all roles associated for each ChEBI ID
        Dictionary of roles ontology, associating for each role, its parent roles
    """
    d_roles_ontology = dict()
    all_roles = dict()
    chebi_ok = 0
    total_nb = len(chebi_ids)
    for chebi_id in chebi_ids:
        roles = set()
        sparql = SPARQLWrapper(endpoint_url)
        sparql.setQuery(f"""
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
        PREFIX owl: <http://www.w3.org/2002/07/owl#>
        PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
        PREFIX dc: <http://purl.org/dc/elements/1.1/>
        PREFIX dcterms: <http://purl.org/dc/terms/>
        PREFIX chebi: <http://purl.obolibrary.org/obo/chebi/>
        PREFIX chebidb: <http://purl.obolibrary.org/obo/CHEBI_>
        PREFIX chebirel: <http://purl.obolibrary.org/obo/chebi#>
        PREFIX obo: <http://purl.obolibrary.org/obo/>
        PREFIX oboInOwl: <http://www.geneontology.org/formats/oboInOwl#>
        PREFIX bp3: <http://www.biopax.org/release/biopax-level3.owl#>
    
        SELECT DISTINCT ?molecule ?moleculeLabel ?role ?roleLabel ?parentRoleLabel
        WHERE {{
            VALUES ?molecule{{ chebidb:{chebi_id }}}                                        
            
            ?molecule rdfs:label ?moleculeLabel .
            
            ?molecule rdfs:subClassOf+ ?restriction .
            ?restriction rdf:type owl:Restriction .
            ?restriction owl:onProperty obo:RO_0000087 .
            ?restriction owl:someValuesFrom/(rdfs:subClassOf*) ?role .
            
            ?role rdfs:subClassOf ?parentRole .
            
            ?parentRole rdfs:label ?parentRoleLabel .
            ?role rdfs:label ?roleLabel .
        }}
        """)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        parent_roles = set()
        for result in results["results"]["bindings"]:
            role = result['roleLabel']['value']
            parent_role = result['parentRoleLabel']['value']
            parent_roles.add(parent_role)
            roles.add(result['roleLabel']['value'])
            if role not in d_roles_ontology:
                d_roles_ontology[role] = set()
            d_roles_ontology[role].add(parent_role)
        if roles:
            chebi_ok += 1
            d_roles_ontology[chebi_id] = list(roles.difference(parent_roles))
            roles.add(ROOTS[CHEBI])
            all_roles[chebi_id] = roles
        else:
            print(f'No ChEBI role found for : {chebi_id}')

    for c, p in d_roles_ontology.items():
        d_roles_ontology[c] = list(p)

    print(f'{chebi_ok}/{total_nb} chebi id with roles associated.')
    return all_roles, d_roles_ontology


def extract_go_classes(go_ids: List[str], endpoint_url: str) \
        -> Tuple[Dict[str, Set[str]], Dict[str, List[str]]]:
    """

    Parameters
    ----------
    go_ids: List[str]
        List of GO IDs
    endpoint_url: str
        URL endpoint string

    Returns
    -------
    Tuple[Dict[str, Set[str]], Dict[str, List[str]]]
        Dictionary of all roles associated for each ChEBI ID
        Dictionary of roles ontology, associating for each role, its parent roles
    """
    d_classes_ontology = dict()
    all_classes = dict()
    for go in go_ids:
        go = go.lower()
        go_classes = set()
        sparql = SPARQLWrapper(endpoint_url)
        sparql.setQuery(f"""
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
        
        SELECT ?goLabel ?parentGoLabel
        WHERE {{
           {go} rdfs:subClassOf* ?go .
           ?go rdfs:label ?goLabel .
           ?go rdf:type owl:Class .
           ?go rdfs:subClassOf ?parentGo .
           ?parentGo rdfs:label ?parentGoLabel .
           ?parentGo rdf:type owl:Class .
        }}
        """)
        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        parent_classes = set()
        for result in results["results"]["bindings"]:
            go_class = result['goLabel']['value']
            parent_class = result['parentGoLabel']['value']
            parent_classes.add(parent_class)
            go_classes.add(go_class)
            go_classes.add(parent_class)
            if parent_class in GO_ROOTS:
                d_classes_ontology[parent_class] = [ROOTS[GO]]
            if go_class not in d_classes_ontology:
                d_classes_ontology[go_class] = set()
            d_classes_ontology[go_class].add(parent_class)
        if go_classes:
            d_classes_ontology[go] = list(go_classes.difference(parent_classes))
            go_classes.add(ROOTS[GO])
            all_classes[go] = go_classes
        else:
            print(f'No GO class found for : {go}')

    for c, p in d_classes_ontology.items():
        d_classes_ontology[c] = list(p)
    return all_classes, d_classes_ontology


def get_parents(child: str, parent_set: Set[str], d_classes_ontology: Dict[str, List[str]],
                root_item) -> Set[str]:
    """ Get recursively from a child class, all its parents classes found in ontology.

    Parameters
    ----------
    child: str
        Child class
    parent_set: Set[str]
        Set of all parents from previous classes
    d_classes_ontology: Dict[str, List[str]]
        Dictionary of the classes ontology of MetaCyc associating for each class its parent classes.
    root_item: str
        Name of the root item of the ontology

    Returns
    -------
    Set[str]
        Set of the union of the set  of child parent classes and the set of all previous parents.
    """
    parents = d_classes_ontology[child]
    for p in parents:
        parent_set.add(p)
    if parents != [root_item]:
        for p in parents:
            parent_set = get_parents(p, parent_set, d_classes_ontology, root_item)
    return parent_set


def get_all_classes(obj_classes: Dict[str, List[str]], d_classes_ontology: Dict[str, List[str]],
                    root_item: str) -> Dict[str, Set[str]]:
    """ Extract all parent classes for each metabolite.

    Parameters
    ----------
    obj_classes: Dict[str, List[str]] (Dict[metabolite, List[class]])
        Dictionary associating for each object the list of +1 parent classes it belongs to.
    d_classes_ontology: Dict[str, List[str]]
        Dictionary of the classes ontology associating for each class its +1 parent classes.
    root_item: str
        Name of the root item of the ontology.

    Returns
    -------
    Dict[str, Set[str]] (Dict[metabolite, Set[class]])
        Dictionary associating for each metabolite the list of all parent classes it belongs to.
    """
    all_classes_met = dict()
    for met, classes in obj_classes.items():
        all_classes = set(classes)
        for c in classes:
            if c != root_item:
                m_classes = get_parents(c, set(d_classes_ontology[c]), d_classes_ontology, root_item)
                all_classes = all_classes.union(m_classes)
        all_classes_met[met] = all_classes
    return all_classes_met


def extract_classes(ontology: str, metabolic_objects: List[str], root: str,
                    d_classes_ontology: Dict[str, List[str]] = None, endpoint_url: str = None)\
        -> Tuple[Dict[str, Set[str]], Dict[str, List[str]]]:
    """ Extract all parent classes (until root) from a list of metabolic objects.

    Parameters
    ----------
    ontology: str
        Name of the ontology considered
    metabolic_objects: List[str]
        List of metabolic object to classify
    root: str
        Root item of the ontology
    d_classes_ontology: Dict[str, List[str]], optional (default=None)
        Dictionary of the classes ontology associating for each class its +1 parent classes.
    endpoint_url: str, optional (default=None)
        URL for the SPARQL server (for GO and ChEBI ontologies)

    Returns
    -------
    Dict[str, Set[str]]
        Dictionary associating to each object all its parent classes (until the root)
    Dict[str, List[str]]
        Dictionary of the classes ontology associating for each class its +1 parent classes.
    """
    if ontology == METACYC or ontology == KEGG or ontology is None:
        leaf_classes = extract_met_classes(metabolic_objects, d_classes_ontology)
        return get_all_classes(leaf_classes, d_classes_ontology, root), d_classes_ontology
    if ontology == EC:
        leaf_classes, d_classes_ontology = extract_ec_classes(metabolic_objects, d_classes_ontology)
        return get_all_classes(leaf_classes, d_classes_ontology, root), d_classes_ontology
    if ontology == CHEBI:
        return extract_chebi_roles(metabolic_objects, endpoint_url)
    if ontology == GO:
        return extract_go_classes(metabolic_objects, endpoint_url)


# For all Ontology - Utils
# --------------------------------------------------------------------------------------------------

def get_abundance_dict(abundances: List[float] or None, metabolic_objects: List[str], ref: bool)\
        -> Dict[str, float]:
    if abundances is None:
        abundances = len(metabolic_objects) * [1]
    if len(metabolic_objects) == len(abundances):
        abundances_dict = {}
        for i in range(len(metabolic_objects)):
            abundances_dict[metabolic_objects[i]] = abundances[i]
    else:
        if ref:
            raise AttributeError(f'Length of "reference_set" parameter must be equal to '
                                 f'"ref_abundances" parameter length : {len(metabolic_objects)} '
                                 f'!= {len(abundances)}')
        else:
            raise AttributeError(f'Length of "metabolic_objects" parameter must be equal to '
                                 f'"abundances" parameter length : {len(metabolic_objects)} '
                                 f'!= {len(abundances)}')
    return abundances_dict


def get_classes_abundance(all_classes: Dict[str, Set[str]], abundances_dict: Dict[str, float],
                          show_leaves: bool) -> Dict[str, float]:
    """ Indicate for each class the number of metabolites found belonging to the class

    Parameters
    ----------
    all_classes: Dict[str, Set[str]] (Dict[metabolite, Set[class]])
        Dictionary associating for each metabolite the list of all parent classes it belongs to.
    abundances_dict: Dict[str, float]
        Dictionary associating for each concept, its abundance value
    show_leaves: bool
        True to show input metabolic objets at sunburst leaves

    Returns
    -------
    Dict[str, float]
        Dictionary associating for each class the weight of concepts found belonging to the class.
    """
    classes_abondance = dict()
    for met, classes in all_classes.items():
        if show_leaves:
            classes_abondance[met] = abundances_dict[met]
        for c in classes:
            if c not in classes_abondance.keys():
                classes_abondance[c] = abundances_dict[met]
            else:
                classes_abondance[c] += abundances_dict[met]
    return dict(reversed(sorted(classes_abondance.items(), key=lambda item: item[1])))


def get_children_dict(parent_dict: Dict[str, List[str]]) -> Dict[str, List[str]]:
    """ Create the children dictionary from the parents dictionary.

    Parameters
    ----------
    parent_dict: Dict[str, List[str]]
        Dictionary associating for each class, its parents classes

    Returns
    -------
    Dict[str, List[str]]
        Dictionary associating for each class, its children classes
    """
    children_dict = dict()
    for c, ps in parent_dict.items():
        for p in ps:
            if p not in children_dict.keys():
                children_dict[p] = list()
            children_dict[p].append(c)
    return children_dict


def reduce_d_ontology(d_classes_ontology, classes_abundance):
    reduced_d_ontology = dict()
    for k, v in d_classes_ontology.items():
        if k in classes_abundance:
            reduced_d_ontology[k] = v
    return reduced_d_ontology
