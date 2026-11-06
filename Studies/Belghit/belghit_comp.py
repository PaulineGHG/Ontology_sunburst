import numpy as np

from ontosunburst.ontosunburst import ontosunburst
import pandas as pd
FILE = 'Belghit2017.csv'
FILE_O = 'Belghit2017-original.csv'


def analysis_1():
    root = 'r'
    with open(FILE, 'r') as f:
        onto = {}
        red_med_d = {}
        red_moy_d = {}
        green_med_d = {}
        green_moy_d = {}
        brown_med_d = {}
        brown_moy_d = {}
        f.__next__()
        f.__next__()
        f.__next__()
        for l in f:
            l = l.split('\t')
            biochem = l[13]
            sub_pwy = l[12]
            sup_pwy = l[11]

            red_med = float(l[14].replace('#VALEUR !', '0').replace(',', '.'))
            green_med = float(l[15].replace('#VALEUR !', '0').replace(',', '.'))
            brown_med = float(l[16].replace('#VALEUR !', '0').replace(',', '.'))

            red_moy = float(l[17].replace('#DIV/0 !', '0').replace(',', '.'))
            green_moy = float(l[18].replace('#DIV/0 !', '0').replace(',', '.'))
            brown_moy = float(l[19].replace('#DIV/0 !', '0').replace(',', '.'))

            onto[biochem] = [sub_pwy]
            onto[sub_pwy] = [sup_pwy]
            onto[sup_pwy] = [root]

            if red_med != 0:
                red_med_d[biochem] = red_med
            if red_moy != 0:
                red_moy_d[biochem] = red_moy
            if green_med != 0:
                green_med_d[biochem] = green_med
            if green_moy != 0:
                green_moy_d[biochem] = green_moy
            if brown_med != 0:
                brown_med_d[biochem] = brown_med
            if brown_moy != 0:
                brown_moy_d[biochem] = brown_moy

    ontosunburst(interest=red_med_d, ontology_dag_input=onto, output='red_med')
    ontosunburst(interest=red_moy_d, ontology_dag_input=onto, output='red_moy')
    ontosunburst(interest=green_med_d, ontology_dag_input=onto, output='green_med')
    ontosunburst(interest=green_moy_d, ontology_dag_input=onto, output='green_moy')
    ontosunburst(interest=brown_med_d, ontology_dag_input=onto, output='brown_med')
    ontosunburst(interest=brown_moy_d, ontology_dag_input=onto, output='brown_moy')


def analysis_data_corrected():
    root = 'r'

    onto = {}

    i_biochem = 13
    i_sub_pw = 12
    i_sup_pw = 11
    i_all = slice(14, 56)
    i_red = slice(14, 28)
    i_green = slice(28, 34)
    i_brown = slice(34, 56)

    red_col = ['red' + str(x+1) for x in range(i_red.stop - i_red.start)]
    green_col = ['green' + str(x + 1) for x in range(i_green.stop - i_green.start)]
    brown_col = ['brown' + str(x + 1) for x in range(i_brown.stop - i_brown.start)]

    df = pd.DataFrame(columns=red_col+green_col+brown_col)
    with open(FILE_O, 'r') as f:
        f.__next__()
        f.__next__()
        f.__next__()
        f.__next__()
        f.__next__()
        f.__next__()
        for l in f:
            l = l[:-1].split('\t')

            # Ontology
            onto[l[i_biochem]] = [l[i_sub_pw]]
            onto[l[i_sub_pw]] = [l[i_sup_pw]]
            onto[l[i_sup_pw]] = [root]

            # DF
            lst = [float("nan") if x == '' else int(x) for x in l[i_all]]
            row = pd.DataFrame([lst], columns=df.columns, index=[l[i_biochem]])
            df = pd.concat([row, df])

    df_n = df.div(df.median(axis=1), axis=0)
    df_n = df_n.astype(float)

    df_n['mean_red'] = df_n[red_col].mean(axis=1)
    df_n['mean_green'] = df_n[green_col].mean(axis=1)
    df_n['mean_brown'] = df_n[brown_col].mean(axis=1)
    df_n['mean_all'] = df_n[red_col + green_col + brown_col].mean(axis=1)

    df_n['med_red'] = df_n[red_col].median(axis=1)
    df_n['med_green'] = df_n[green_col].median(axis=1)
    df_n['med_brown'] = df_n[brown_col].median(axis=1)
    df_n['med_all'] = df_n[red_col + green_col + brown_col].median(axis=1)

    red_med_d = df_n['med_red'].dropna().to_dict()
    red_moy_d = df_n['mean_red'].dropna().to_dict()
    green_med_d = df_n['med_green'].dropna().to_dict()
    green_moy_d = df_n['mean_green'].dropna().to_dict()
    brown_med_d = df_n['med_brown'].dropna().to_dict()
    brown_moy_d = df_n['mean_brown'].dropna().to_dict()
    all_med_d = df_n['med_all'].dropna().to_dict()
    all_moy_d = df_n['mean_all'].dropna().to_dict()

    ontosunburst(interest=red_med_d, ontology_dag_input=onto, output='red_med')
    ontosunburst(interest=red_moy_d, ontology_dag_input=onto, output='red_moy')
    ontosunburst(interest=green_med_d, ontology_dag_input=onto, output='green_med')
    ontosunburst(interest=green_moy_d, ontology_dag_input=onto, output='green_moy')
    ontosunburst(interest=brown_med_d, ontology_dag_input=onto, output='brown_med')
    ontosunburst(interest=brown_moy_d, ontology_dag_input=onto, output='brown_moy')

    mul = 10000
    ontosunburst(interest={k: int(v*mul) for k, v in red_med_d.items()},
                 reference={k: int(v*mul) for k, v in all_med_d.items()},
                 colorization='e', ontology_dag_input=onto, output='red_med_e')
    ontosunburst(interest={k: int(v*mul) for k, v in red_moy_d.items()},
                 reference={k: int(v*mul) for k, v in all_moy_d.items()},
                 colorization='e', ontology_dag_input=onto, output='red_moy_e')
    ontosunburst(interest={k: int(v*mul) for k, v in green_med_d.items()},
                 reference={k: int(v*mul) for k, v in all_med_d.items()},
                 colorization='e', ontology_dag_input=onto, output='green_med_e')
    ontosunburst(interest={k: int(v*mul) for k, v in green_moy_d.items()},
                 reference={k: int(v*mul) for k, v in all_moy_d.items()},
                 colorization='e', ontology_dag_input=onto, output='green_moy_e')
    ontosunburst(interest={k: int(v*mul) for k, v in brown_med_d.items()},
                 reference={k: int(v*mul) for k, v in all_med_d.items()}, 
                 colorization='e', ontology_dag_input=onto, output='brown_med_e')
    ontosunburst(interest={k: int(v*mul) for k, v in brown_moy_d.items()},
                 reference={k: int(v*mul) for k, v in all_moy_d.items()},
                 colorization='e', ontology_dag_input=onto, output='brown_moy_e')

    ontosunburst(interest=red_med_d, reference=all_med_d,
                 colorization='d', ontology_dag_input=onto, output='red_med_d')
    ontosunburst(interest=red_moy_d, reference=all_moy_d,
                 colorization='d', ontology_dag_input=onto, output='red_moy_d')
    ontosunburst(interest=green_med_d, reference=all_med_d,
                 colorization='d', ontology_dag_input=onto, output='green_med_d')
    ontosunburst(interest=green_moy_d, reference=all_moy_d,
                 colorization='d', ontology_dag_input=onto, output='green_moy_d')
    ontosunburst(interest=brown_med_d, reference=all_med_d,
                 colorization='d', ontology_dag_input=onto, output='brown_med_d')
    ontosunburst(interest=brown_moy_d, reference=all_moy_d,
                 colorization='d', ontology_dag_input=onto, output='brown_moy_d')

    # df_log2 = np.log2(df_normalized)


analysis_data_corrected()
