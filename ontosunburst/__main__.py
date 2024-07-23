from ontosunburst.ontosunburst import *
import argparse


def get_command_line_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, required=True, help='Interest set input')
    parser.add_argument('--ref', type=str, required=False, help='Reference set input')
    parser.add_argument('--ontology', type=str, required=False, help='Ontology used')
    parser.add_argument('--root', type=str, required=False, help='Ontology root')
    parser.add_argument('--analysis', type=str, required=False, help='Type of analysis')
    parser.add_argument('--output', type=str, required=False, help='Output path+name')
    parser.add_argument('--class_ontology', type=str, required=False,
                        help='Class ontology json file')
    parser.add_argument('--labels', type=str, required=False, help='Labels json file')
    parser.add_argument('--url', type=str, required=False, help='Endpoint URL (SPARQL server)')
    parser.add_argument('--test', type=str, required=False, help='Enrichment stat test')
    parser.add_argument('--total', action='store_true', required=False, help='Total branch value')
    parser.add_argument('--rcut', type=str, required=False, help='Type of root cut')
    parser.add_argument('--r_base', action='store_true', required=False, help='Reference base')
    parser.add_argument('--show_leaves', action='store_true', required=False, help='Show leaves')
    parser.add_argument('--total', action='store_true', required=False, help='Total branch value')
    args = parser.parse_args()
    return args


def main():
    args = get_command_line_args()
    metabolic_objects, abundances = extract_input(args.input)
    reference_set, ref_abundances = extract_input(args.ref)
    ontosunburst(metabolic_objects=metabolic_objects,
                 ontology=args.ontology,
                 root=args.root,
                 abundances=abundances,
                 reference_set=reference_set,
                 ref_abundances=ref_abundances,
                 analysis=args.analysis,
                 output=args.output,
                 write_output=True,
                 class_ontology=args.class_ontology,
                 labels=args.labels,
                 endpoint_url=args.url,
                 test=args.test,
                 total=args.total,
                 root_cut=args.rcut,
                 ref_base=args.r_base,
                 show_leaves=args.show_leaves)


def extract_input(input_file):
    id_lst = []
    ab_lst = []
    with open(input_file, 'r') as f:
        for l in f:
            l = l.strip().split('\n')
            id_lst.append(l[0])
            if len(l) == 2:
                ab_lst.append(l[1])
            else:
                ab_lst.append(1)
    return id_lst, ab_lst
