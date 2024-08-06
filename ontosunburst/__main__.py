from ontosunburst.ontosunburst import *
import argparse


def get_command_line_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', '-i', type=str, required=True, help='Interest set input')
    parser.add_argument('--ref', '-ir', type=str, required=False, help='Reference set input')
    parser.add_argument('--ontology', '--onto', type=str, required=False, help='Ontology used')
    parser.add_argument('--root', '-r', type=str, required=False, help='Ontology root')
    parser.add_argument('--analysis', '-a', type=str, required=False, help='Type of analysis')
    parser.add_argument('--output', '-o', type=str, required=False, help='Output path+name')
    parser.add_argument('--class_ontology', '-cl', type=str, required=False,
                        help='Class ontology json file')
    parser.add_argument('--labels', '-l', type=str, required=False, help='Labels json file')
    parser.add_argument('--url', type=str, required=False, help='Endpoint URL (SPARQL server)')
    parser.add_argument('--test', '-t', type=str, required=False, help='Enrichment stat test')
    parser.add_argument('--rcut', type=str, required=False, help='Type of root cut')
    parser.add_argument('--pcut', type=str, required=False, help='Type of path cut')
    parser.add_argument('--r_base', action='store_true', required=False, help='Reference base')
    parser.add_argument('--show_leaves', '-sl',  action='store_true', required=False,
                        help='Show leaves')
    parser.add_argument('--kwargs', nargs=argparse.REMAINDER, help="Additional keyword arguments")
    args = parser.parse_args()
    return args


def main():
    args = get_command_line_args()
    kwargs = {}
    if args.kwargs:
        for arg in args.kwargs:
            if "=" in arg:
                key, value = arg.split("=", 1)
                kwargs[key] = value
            else:
                raise ValueError(f"Argument {arg} is not in the form key=value")

    metabolic_objects, abundances = extract_input(args.input)
    reference_set, ref_abundances = extract_input(args.ref)
    ontosunburst(interest_set=metabolic_objects,
                 ontology=args.ontology,
                 root=args.root,
                 abundances=abundances,
                 scores={},
                 reference_set=reference_set,
                 ref_abundances=ref_abundances,
                 analysis=args.analysis,
                 output=args.output,
                 write_output=True,
                 class_ontology=args.class_ontology,
                 labels=args.labels,
                 endpoint_url=args.url,
                 test=args.test,
                 root_cut=args.rcut,
                 path_cut='',
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
