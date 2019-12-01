# import pandas as pd
#
# a = pd.read_csv('renta_demografía_2016/renta_media.csv', skiprows=4, sep=',', error_bad_lines=False)
# b = pd.read_csv('renta_demografía_2016/indicadores_demográficos.csv', skiprows=5, sep=',', error_bad_lines=False)
#
#
# print(a)

import json

archivo = "D:/Descargas/SECC_CPV_E_20111101_01_R_INE.json"
with open(archivo, encoding='utf-8') as json_file:
    data = json.load(json_file)

counter = 0
for i in data:
    if counter < 10:
        print('HOLAAAAA   ', i)
        counter += 1
    else:
        break