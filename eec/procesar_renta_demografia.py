"""
En este código vamos a crear las tablas de cada uno de
los archivos csv de renta que tenemos. El archivo incluye información de 2015 y 2016. La de 2015 la desechamos,
ordenamos las columnas, y reajustamos los municipios. Por último, llenamos valores de casillas que estén vacías con la
información de la región, y si no existe, se queda en NaN.
"""

import pandas as pd
import numpy as np
from tqdm import tqdm


def procesar_renta(carpeta, años):
    """
    Crea un dataframe con toda la información de los datos de la renta. En el caso de 2015 2016, hay una separación
    de cada una de las columnas en dos años. La variable años sirve para indicar que esas dos columnas existen y se
    va a generar un MultiIndex basado en esa composición anual.
    :param carpeta:
    :param años:
    :return:
    """
    df = pd.read_csv(carpeta + '/renta_media.csv', skiprows=4, sep=';', error_bad_lines=False, dtype=str)
    df_name = pd.DataFrame(data=df.iloc[1:-4,0].values, columns=['nombre'])
    df_rest = df.iloc[1:-4,1:-1] # La primera columna es el codigo, y la última es NaNs. Las filas eliminadas son NaNs también
    
    df_rest.columns = pd.MultiIndex.from_product([['Renta media por persona', 'Renta media por hogar', 
                                                   'Mediana de renta por unidad de consumo', 'Renta bruta media por persona', 
                                                   'Renta bruta media por hogar'], años],
                                                 names=['renta', 'año'])
    
    codigo = [n.split(' ')[0] for n in df_name['nombre'].values]
    
    """
    Ahora que las columnas están completas, vamos a sustituir las celdas con puntos con np.NaN, por practicidad.
    Antes de eso, las celdas que originariamente tienen nan son porque no hay información de la sección, pero sí del municipio
    luego a esas celdas se les añade la información del municipio. Las celdas con punto no tienen información ni del municipio,
    luego se les añade el np.NaN.
    """
    
    df_name.index, df_rest.index = codigo, codigo
    df_rest[df_rest == '.'] = 0
    df_rest = df_rest.astype(float)

    nanpos = df_rest[pd.isnull(df_rest).all(axis=1)].index # Si se pone np.any se va gipuzkoa porque en 2015 no tenia renta

    for row in tqdm(nanpos):
        if len(row) < 6:
            df_rest.loc[row, :] = np.NaN
        else:
            df_rest.loc[row, :] = df_rest.loc[row[:5], :]

    df_rest[df_rest == 0] = np.NaN
    df_rest.to_pickle(carpeta + '/renta_media_procesado.pickle')


def procesar_demografia(fichero, columnas, cifras_cutoff=8):
    """
    
    """
    df = pd.read_csv(fichero, skiprows=5, sep=';', error_bad_lines=False, names=None, header=None, dtype=str)
    
    df_name = pd.DataFrame(data=df.iloc[len(columnas):-4, 0].values, columns=['nombre'])
    df_rest = df.iloc[len(columnas):-4, 1:-1]

    
    df_rest.columns = pd.MultiIndex.from_product([list(dict.fromkeys([x for x in df.loc[row, :].values.tolist()
                                                                      if isinstance(x, str)])) for row in
                                                  range(len(columnas))], names=columnas)
    codigo = [n.split(' ')[0] for n in df_name['nombre'].values]

    df_name.index, df_rest.index = codigo, codigo
    df_rest[df_rest == '.'] = np.NaN
    df_rest[df_rest == '..'] = np.NaN
    # En estos datos no hay información de sección, pero sí de distrito. Así, primero quitamos las secciones y luego rellenamos
    # la información de los distritos. CUIDADO! Para el archivo de datos demográficos generales SI hay info de todo. Así que ese nos lo tenemos que quedar.
    df_rest = df_rest.loc[[i for i in codigo if len(i) < cifras_cutoff]].apply(pd.to_numeric, errors='coerce')
    nan_rows = df_rest.loc[np.all(np.isnan(df_rest), 1)].index

    for row in tqdm(nan_rows):
        df_rest.loc[row, :] = df_rest.loc[row[:5], :]

    df_rest.to_pickle(fichero.replace('.csv', '.pickle'))
    
    
def reemplazar_puntos_comas_csv(archivo):
    with open(archivo, 'r') as fi:
        texto = fi.read()
        texto = texto.replace('.', '').replace(',', '.')
    
    with open(archivo, 'w') as fo:
        fo.write(texto)