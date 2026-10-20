from ontosunburst.ontosunburst import ontosunburst
FILE = 'Belghit2017.csv'

root = 'r'
with open(FILE, 'r') as f:
    onto = {}
    red_med_d = {}
    red_moy_d = {}
    green_med_d = {}
    green_moy_d = {}
    brown_med_d = {}
    brown_moy_d = {}
    all_med_d = {}
    all_moy_d = {}
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

        all_med = float(l[20].replace('#VALUE!', '0').replace(',', '.'))
        all_moy = float(l[21].replace('#DIV/0!', '0').replace(',', '.'))

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
        if all_med != 0:
            all_med_d[biochem] = all_med
        if all_moy != 0:
            all_moy_d[biochem] = all_moy

# ontosunburst(interest=red_med_d, ontology_dag_input=onto, output='red_med')
# ontosunburst(interest=red_moy_d, ontology_dag_input=onto, output='red_moy')
# ontosunburst(interest=green_med_d, ontology_dag_input=onto, output='green_med')
# ontosunburst(interest=green_moy_d, ontology_dag_input=onto, output='green_moy')
# ontosunburst(interest=brown_med_d, ontology_dag_input=onto, output='brown_med')
# ontosunburst(interest=brown_moy_d, ontology_dag_input=onto, output='brown_moy')

# ontosunburst(interest=red_med_d, reference=all_med_d, ontology_dag_input=onto, output='red_med_e', colorization='e')
# ontosunburst(interest=red_moy_d, reference=all_moy_d, ontology_dag_input=onto, output='red_moy_e', colorization='e')
# ontosunburst(interest=green_med_d, reference=all_med_d, ontology_dag_input=onto, output='green_med_e', colorization='e')
# ontosunburst(interest=green_moy_d, reference=all_moy_d, ontology_dag_input=onto, output='green_moy_e', colorization='e')
# ontosunburst(interest=brown_med_d, reference=all_med_d, ontology_dag_input=onto, output='brown_med_e', colorization='e')
# ontosunburst(interest=brown_moy_d, reference=all_moy_d, ontology_dag_input=onto, output='brown_moy_e', colorization='e')

ontosunburst(interest=red_med_d, reference=all_med_d, ontology_dag_input=onto, output='red_med_d', colorization='d')
ontosunburst(interest=red_moy_d, reference=all_moy_d, ontology_dag_input=onto, output='red_moy_d', colorization='d')
ontosunburst(interest=green_med_d, reference=all_med_d, ontology_dag_input=onto, output='green_med_d', colorization='d')
ontosunburst(interest=green_moy_d, reference=all_moy_d, ontology_dag_input=onto, output='green_moy_d', colorization='d')
ontosunburst(interest=brown_med_d, reference=all_med_d, ontology_dag_input=onto, output='brown_med_d', colorization='d')
ontosunburst(interest=brown_moy_d, reference=all_moy_d, ontology_dag_input=onto, output='brown_moy_d', colorization='d')