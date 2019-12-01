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
    df = pd.read_csv(carpeta + '/renta_media.csv', skiprows=4, sep=';', error_bad_lines=False)

    df_name = pd.DataFrame(data=df.iloc[1:-4,0].values, columns=['nombre'])
    df_rest = df.iloc[1:-4,1:-1]

    df_rest.columns = pd.MultiIndex.from_product([['Renta media por persona', 'Renta media por hogar'], años],
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

    nanpos = df_rest[pd.isnull(df_rest).any(axis=1)].index

    for row in tqdm(nanpos):
        if len(row) < 6:
            df_rest.loc[row, :] = np.NaN
        else:
            df_rest.loc[row, :] = df_rest.loc[row[:5], :]

    df_rest[df_rest == 0] = np.NaN
    df_rest.to_pickle(carpeta + '/renta_media_procesado.pickle')


def procesar_demografia(fichero, columnas):
    df = pd.read_csv(fichero, skiprows=4, sep=';', error_bad_lines=False, names=None, header=None)

    df_name = pd.DataFrame(data=df.iloc[len(columnas):-4, 0].values, columns=['nombre'])
    df_rest = df.iloc[len(columnas):-4, 1:-1]

    df_rest.columns = pd.MultiIndex.from_product([list(dict.fromkeys([x for x in df.loc[row, :].values.tolist()
                                                                      if isinstance(x, str)])) for row in
                                                  range(len(columnas))],
                                                 names=columnas)

    codigo = [n.split(' ')[0] for n in df_name['nombre'].values]

    df_name.index, df_rest.index = codigo, codigo
    df_rest[df_rest == '.'] = np.NaN
    df_rest[df_rest == '..'] = np.NaN

    # En estos datos no hay información de sección, pero sí de distrito. Así, primero quitamos las secciones y luego rellenamos
    # la información de los distritos
    df_rest = df_rest.loc[[i for i in codigo if len(i) < 8]].astype(float)

    nan_rows = df_rest.loc[np.all(np.isnan(df_rest), 1)].index

    for row in tqdm(nan_rows):
        df_rest.loc[row, :] = df_rest.loc[row[:5], :]

    df_rest.to_pickle(fichero.replace('.csv', '.pickle'))