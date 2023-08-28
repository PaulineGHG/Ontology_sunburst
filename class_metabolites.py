mport json

from padmet.classes import PadmetRef
from typing import Set

from obj_extraction import *
from ontology import *
from sunburst_fig import *


# WORKFLOW ==========================================================================================================

def proportion_workflow(metabolites: Set[str], class_file: str, padmet_ref: str, output: str, full=True):
    """ Classify and plot a list of metabolites (MetaCyc Ids)

    Parameters
    ----------
    metabolites: Set[str]
        Set of metabolites to classify
    class_file: str
        Path to class json file
    padmet_ref: str
        Path to reference padmet file from MetaCyc.
    output: str
        Path to output to save figure.
    full: bool
    """
    with open(class_file, 'r') as f:
        d_classes_ontology = json.load(f)
    padmet_ref = PadmetRef(padmet_ref)
    # d_classes_ontology = get_dict_ontology(padmet_ref)
    met_classes = extract_classes(metabolites, padmet_ref)
    all_classes = get_all_classes(met_classes, d_classes_ontology, 'FRAMES')
    write_met_classes(all_classes, output)
    classes_abundance = get_classes_abondance(all_classes)
    data = get_fig_parameters(classes_abundance, d_classes_ontology,
                              get_children_dict(d_classes_ontology), 'FRAMES', full)
    data = get_data_proportion(data)
    generate_sunburst_fig(data, output, 'proportion')


def comparison_workflow(metabolites_interest, metabolites_base, class_file, padmet_ref, output,  test='Binomial',
                        full=True):
    with open(class_file, 'r') as f:
        d_classes_ontology = json.load(f)
    padmet_ref = PadmetRef(padmet_ref)
    i_met_classes = extract_classes(metabolites_interest, padmet_ref)
    i_all_classes = get_all_classes(i_met_classes, d_classes_ontology, 'FRAMES')
    i_classes_abondance = get_classes_abondance(i_all_classes)
    b_met_classes = extract_classes(metabolites_base, padmet_ref)
    b_all_classes = get_all_classes(b_met_classes, d_classes_ontology, 'FRAMES')
    b_classes_abundance = get_classes_abondance(b_all_classes)
    data = get_fig_parameters(i_classes_abondance, d_classes_ontology,
                              get_children_dict(d_classes_ontology), 'FRAMES', full)
    data = get_data_prop_diff(data, b_classes_abundance)
    generate_sunburst_fig(data, output, 'comparison', b_classes_abundance, test)


def all_workflow(metabolites_interest, metabolites_base, class_file, padmet_ref, output, full=True):
    with open(class_file, 'r') as f:
        d_classes_ontology = json.load(f)
    padmet_ref = PadmetRef(padmet_ref)
    i_met_classes = extract_classes(metabolites_interest, padmet_ref)
    i_all_classes = get_all_classes(i_met_classes, d_classes_ontology, 'FRAMES')
    i_classes_abondance = get_classes_abondance(i_all_classes)
    b_met_classes = extract_classes(metabolites_base, padmet_ref)
    b_all_classes = get_all_classes(b_met_classes, d_classes_ontology, 'FRAMES')
    b_classes_abundance = get_classes_abondance(b_all_classes)
    data = get_fig_parameters(i_classes_abondance, d_classes_ontology,
                              get_children_dict(d_classes_ontology), 'FRAMES', full)
    data = get_data_proportion(data)
    data = get_data_prop_diff(data, b_classes_abundance)
    generate_dash_interactive_fig(data, output)


def pathways_workflow_proportion(pw_classes, class_file, output, full=True):
    with open(class_file, 'r') as f:
        d_classes_ontology = json.load(f)
    # d_classes_ontology = get_dict_ontology(padmet_ref)
    all_classes = get_all_classes(pw_classes, d_classes_ontology, 'FRAMES')
    # write_met_classes(all_classes, output)
    classes_abondance = get_classes_abondance(all_classes)
    data = get_fig_parameters(classes_abondance, d_classes_ontology,
                              get_children_dict(d_classes_ontology), 'FRAMES', full)
    data = get_data_proportion(data)
    generate_sunburst_fig(data, output, 'proportion')


def pathways_workflow_comparison(pw_cls_interest, pw_cls_base, class_file, output,  test='Binomial',
                        full=True):
    with open(class_file, 'r') as f:
        d_classes_ontology = json.load(f)
    i_all_classes = get_all_classes(pw_cls_interest, d_classes_ontology, 'FRAMES')
    i_classes_abondance = get_classes_abondance(i_all_classes)
    b_all_classes = get_all_classes(pw_cls_base, d_classes_ontology, 'FRAMES')
    b_classes_abundance = get_classes_abondance(b_all_classes)
    data = get_fig_parameters(i_classes_abondance, d_classes_ontology,
                              get_children_dict(d_classes_ontology), 'FRAMES', full)
    data = get_data_prop_diff(data, b_classes_abundance)
    generate_sunburst_fig(data, output, 'comparison', b_classes_abundance, test)


def chebi_roles_workflow_proportion(chebi_ids, endpoint_url, output, full=True):
    all_classes, d_roles_ontology = extract_chebi_roles(chebi_ids, endpoint_url)
    classes_abondance = get_classes_abondance(all_classes)
    data = get_fig_parameters(classes_abondance, d_roles_ontology,
                              get_children_dict(d_roles_ontology), 'role', full)
    data = get_data_proportion(data)
    generate_sunburst_fig(data, output, 'proportion')


def chebi_roles_workflow_comparison(chebi_interest, chebi_base, endpoint_url, output,  test='Binomial',
                                    full=True):
    i_all_roles, d_roles_ontology = extract_chebi_roles(chebi_interest, endpoint_url)
    i_roles_abondance = get_classes_abondance(i_all_roles)
    b_all_roles, d_roles_ontology = extract_chebi_roles(chebi_base, endpoint_url)
    b_roles_abundance = get_classes_abondance(b_all_roles)
    data = get_fig_parameters(i_roles_abondance, d_roles_ontology,
                              get_children_dict(d_roles_ontology), 'FRAMES', full)
    data = get_data_prop_diff(data, b_roles_abundance)
    generate_sunburst_fig(data, output, 'comparison', b_roles_abundance, test)


def go_workflow_proportion(go_abundance, class_file, output, full=True):
    with open(class_file, 'r') as f:
        d_classes_ontology = json.load(f)
    for go in go_abundance:
        if go not in d_classes_ontology:
            print(go)
    data = get_fig_parameters(go_abundance, d_classes_ontology,
                              get_children_dict(d_classes_ontology), 'FRAMES', full)
    data = get_data_proportion(data)
    generate_sunburst_fig(data, output, 'proportion')


def ec_workflow_proportion(ec_classes, class_file, names_file, output, full=True):
    with open(class_file, 'r') as f:
        d_classes_ontology = json.load(f)
    with open(names_file, 'r') as f:
        names = json.load(f)
    all_classes = get_all_classes(ec_classes, d_classes_ontology, 'Enzyme')
    classes_abondance = get_classes_abondance(all_classes)
    data = get_fig_parameters(classes_abondance, d_classes_ontology,
                              get_children_dict(d_classes_ontology), 'Enzyme', full, names)
    data = get_data_proportion(data)
    generate_sunburst_fig(data, output, 'proportion')


# EXTRAS ==============================================================================================================

def write_met_classes(all_classes, output):
    with open(f'{output}.tsv', 'w') as f:
        for met, classes, in all_classes.items():
            f.write(f'{met}\t{classes}\n')

