"""
En este código vamos a aunar todas las provincias, y vamos a eliminar la información de municipios y distritos, quedándonos
sólo c
"""

import pandas as pd
import numpy as np
from tqdm import tqdm
import os


def procesar_datos_2011(directorio):
    lista_archivos_ccaa = [i for i in os.listdir(directorio) if '_ccaa' in i]

    lista_df_ccaas = []

    for ccaa in lista_archivos_ccaa:
        df = pd.read_csv(directorio + '/' + ccaa, sep=',')
        lista_df_ccaas.append(df)

    df_indicadores = pd.read_csv(directorio + '/' + 'indicadores.csv', sep=';', header=None, encoding='latin1',
                                 error_bad_lines=False)
    dict_indicadores = dict(zip(df_indicadores.iloc[:, 0], df_indicadores.iloc[:, 1]))
    df_ccaas_entero = pd.concat(lista_df_ccaas)
    df_ccaas_entero.columns = [dict_indicadores[i] if i in dict_indicadores.keys() else i for i in df_ccaas_entero.columns.values]

    df_ccaas_entero.to_pickle(directorio + '/' + 'censo_2011.pickle')


