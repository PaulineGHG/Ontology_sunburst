import json
import os.path
from padmet.classes.padmetRef import PadmetRef

CLASSES_SUFFIX = '_classes.json'
NAMES_SUFFIX = '_names.json'
EC_ROOT = 'Enzyme'


# METACYC
# ==================================================================================================
def generate_metacyc_input(input_padmet, version, output_name='MetaCyc'):
    output = output_name + version + CLASSES_SUFFIX
    pref = PadmetRef(input_padmet)
    rels = pref.getAllRelation()
    classes = dict()
    for r in rels:
        if r.type == 'is_a_class':
            if r.id_in not in classes:
                classes[r.id_in] = set()
            classes[r.id_in].add(r.id_out)
    for c, p in classes.items():
        classes[c] = list(p)
    with open(output, 'w') as o:
        json.dump(fp=o, obj=classes, indent=1)


# EC
# ==================================================================================================
def generate_ec_input(enzclass_txt, enzyme_dat, version, output='enzyme'):
    output_name = output + version + NAMES_SUFFIX
    output_class = output + version + CLASSES_SUFFIX
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
        return EC_ROOT
    else:
        ec_lst = ec.split('.')
        ec_lst[-(ec_lvl+1)] = '-'
        return '.'.join(ec_lst)

# KEGG
# ==================================================================================================
