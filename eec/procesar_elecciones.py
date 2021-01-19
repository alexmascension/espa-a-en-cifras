"""
En este código vamos a crear las tablas de cada uno de
los archivos DAT en un formato más legible para nosotros,
y luego poder cruzarlas entre sí.
La información de cada fichero está en el .doc, que es
el que parsearemos de acuerdo con la información que
nos proporcione.
"""

"""
Por norma
- Los nombres de las columnas irán en minúscula
"""
import pandas as pd
from tqdm import tqdm
import numpy as np

carpeta = 'congreso_2019_04/'
sufijo = '021904'

def procesar_elecciones(carpeta, sufijo):
    print('Procesando archivo 3')
    # 3. CANDIDATURAS
    dict_3 = {'eleccion': [],
              'año': [],
              'mes': [],
              'codigo candidatura': [],
              'siglas candidatura': [],
              'denominacion candidatura': [],
              'codigo provincial': [],
              'codigo autonomico': [],
              'codigo nacional': []}

    with open(carpeta + '03' + sufijo + '.DAT', 'r', encoding="latin-1") as f:
        for l in f.readlines():
            dict_3['eleccion'].append(l[0:2])
            dict_3['año'].append(l[2:6])
            dict_3['mes'].append(l[6:8])
            dict_3['codigo candidatura'].append(l[8:14])
            dict_3['siglas candidatura'].append(l[14:64])
            dict_3['denominacion candidatura'].append(l[64:214])
            dict_3['codigo provincial'].append(l[214:220])
            dict_3['codigo autonomico'].append(l[220:226])
            dict_3['codigo nacional'].append(l[226:232])

    df_3 = pd.DataFrame(dict_3)
    df_3.to_csv(carpeta + '03' + sufijo + '.csv', index=None, sep=';')

    # 4. CANDIDATOS
    print('Procesando archivo 4')
    dict_4 = {'eleccion': [],
              'año': [],
              'mes': [],
              'numero de vuelta': [],
              'codigo provincia': [],
              'districto electoral': [],
              'codigo municipio': [],
              'codigo codigo candidatura': [],
              'orden candidato': [],
              'tipo candidato': [],
              'nombre': [],
              'primer apellido': [],
              'segundo apellido': [],
              'sexo': [],
              'dia nacimiento': [],
              'mes nacimiento': [],
              'año nacimiento': [],
              'dni': [],
              'candidato elegido': [], }

    with open(carpeta + '04' + sufijo + '.DAT', 'r', encoding="latin-1") as f:
        for l in f.readlines():
            dict_4['eleccion'].append(l[:2])
            dict_4['año'].append(l[2:6])
            dict_4['mes'].append(l[6:8])
            dict_4['numero de vuelta'].append(l[8:9])
            dict_4['codigo provincia'].append(l[9:11])
            dict_4['districto electoral'].append(l[11:12])
            dict_4['codigo municipio'].append(l[12:15])
            dict_4['codigo codigo candidatura'].append(l[15:21])
            dict_4['orden candidato'].append(l[21:24])
            dict_4['tipo candidato'].append(l[24:25])
            dict_4['nombre'].append(l[25:50])
            dict_4['primer apellido'].append(l[50:75])
            dict_4['segundo apellido'].append(l[75:100])
            dict_4['sexo'].append(l[100:101])
            dict_4['dia nacimiento'].append(l[101:103])
            dict_4['mes nacimiento'].append(l[103:105])
            dict_4['año nacimiento'].append(l[105:109])
            dict_4['dni'].append(l[109:119])
            dict_4['candidato elegido'].append(l[119:120])

    df_4 = pd.DataFrame(dict_4)
    df_4.to_csv(carpeta + '04' + sufijo + '.csv', index=None, sep=';')

    # 5.MUNICIPIOS
    print('Procesando archivo 5')
    dict_5 = {'tipo eleccion': [], 'año eleccion': [], 'mes eleccion': [], 'numero vuelta': [],
              'codigo comunidad autonoma': [],
              'codigo provincia': [], 'codigo municipio': [], 'numero distrito': [], 'nombre municipio': [],
              'codigo distrito': [],
              'codigo partido judicial': [], 'codigo diputacion provincial': [], 'codigo comarca': [], 'poblacion': [],
              'numero mesas': [],
              'censo ine': [], 'censo escrutinio': [], 'censo cere': [], 'votantes cere': [],
              'votances primer avance': [], 'votantes segundo avance': [],
              'votos blancos': [], 'votos nulos': [], 'votos candidaturas': [], 'numero escaños': [],
              'votos afirmativos referendum': [],
              'votos negativos referendum': [], 'datos oficiales': []}
    dict_5_types = [int] * 8 + [str] + [int] * 18 + [str]

    with open(carpeta + '05' + sufijo + '.DAT', 'r', encoding="latin-1") as f:
        for l in f.readlines():
            dict_5['tipo eleccion'].append(l[:2])
            dict_5['año eleccion'].append(l[2:6])
            dict_5['mes eleccion'].append(l[6:8])
            dict_5['numero vuelta'].append(l[8:9])
            dict_5['codigo comunidad autonoma'].append(l[9:11])
            dict_5['codigo provincia'].append(l[11:13])
            dict_5['codigo municipio'].append(l[13:16])
            dict_5['numero distrito'].append(l[16:18])
            dict_5['nombre municipio'].append(l[18:118].strip())
            dict_5['codigo distrito'].append(l[118:119])
            dict_5['codigo partido judicial'].append(l[119:122])
            dict_5['codigo diputacion provincial'].append(l[122:125])
            dict_5['codigo comarca'].append(l[125:128])
            dict_5['poblacion'].append(l[128:136])
            dict_5['numero mesas'].append(l[136:141])
            dict_5['censo ine'].append(l[141:149])
            dict_5['censo escrutinio'].append(l[149:157])
            dict_5['censo cere'].append(l[157:165])
            dict_5['votantes cere'].append(l[165:173])
            dict_5['votances primer avance'].append(l[173:181])
            dict_5['votantes segundo avance'].append(l[181:189])
            dict_5['votos blancos'].append(l[189:197])
            dict_5['votos nulos'].append(l[197:205])
            dict_5['votos candidaturas'].append(l[205:213])
            dict_5['numero escaños'].append(l[213:216])
            dict_5['votos afirmativos referendum'].append(l[216:224])
            dict_5['votos negativos referendum'].append(l[224:232])
            dict_5['datos oficiales'].append(l[232:233])

    df_5 = pd.DataFrame(dict_5)
    for col in range(len(dict_5_types)):
        col_name = df_5.columns[col]
        df_5[col_name] = df_5[col_name].astype(dict_5_types[col])

    df_5.to_csv(carpeta + '05' + sufijo + '.csv', index=None, sep=';')

    # 6. CANDIDATURAS MUNICIPIOS
    print('Procesando archivo 6')
    dict_6 = {'tipo eleccion': [], 'año eleccion': [], 'mes eleccion': [], 'numero vuelta': [],
              'codigo provincia': [], 'codigo municipio': [], 'numero distrito': [], 'codigo candidatura': [],
              'votos candidatura': [], 'numero candidatos': []}
    dict_6_types = [int] * 10

    with open(carpeta + '06' + sufijo + '.DAT', 'r', encoding="latin-1") as f:
        for l in f.readlines():
            dict_6['tipo eleccion'].append(l[:2])
            dict_6['año eleccion'].append(l[2:6])
            dict_6['mes eleccion'].append(l[6:8])
            dict_6['numero vuelta'].append(l[8:9])
            dict_6['codigo provincia'].append(l[9:11])
            dict_6['codigo municipio'].append(l[11:14])
            dict_6['numero distrito'].append(l[14:16])
            dict_6['codigo candidatura'].append(l[16:22])
            dict_6['votos candidatura'].append(l[22:30])
            dict_6['numero candidatos'].append(l[30:33])

    df_6 = pd.DataFrame(dict_6)
    for col in range(len(dict_6_types)):
        col_name = df_6.columns[col]
        df_6[col_name] = df_6[col_name].astype(dict_6_types[col])
    df_6.to_csv(carpeta + '06' + sufijo + '.csv', index=None, sep=';')

    # 7. DATOS MUNICIPIO
    print('Procesando archivo 7')
    dict_7 = {'tipo eleccion': [], 'año eleccion': [], 'mes eleccion': [], 'numero vuelta': [],
              'codigo comunidad autonoma': [],
              'codigo provincia': [], 'codigo distrito': [], 'nombre territorio': [], 'poblacion': [],
              'numero mesas': [],
              'censo ine': [], 'censo escrutinio': [], 'censo cere': [], 'votantes cere': [],
              'votances primer avance': [], 'votantes segundo avance': [],
              'votos blancos': [], 'votos nulos': [], 'votos candidaturas': [], 'numero escaños': [],
              'votos afirmativos referendum': [],
              'votos negativos referendum': [], 'datos oficiales': []}
    dict_7_types = [int] * 7 + [str] + [int] * 14 + [str]

    with open(carpeta + '07' + sufijo + '.DAT', 'r', encoding="latin-1") as f:
        for l in f.readlines():
            dict_7['tipo eleccion'].append(l[:2])
            dict_7['año eleccion'].append(l[2:6])
            dict_7['mes eleccion'].append(l[6:8])
            dict_7['numero vuelta'].append(l[8:9])
            dict_7['codigo comunidad autonoma'].append(l[9:11])
            dict_7['codigo provincia'].append(l[11:13])
            dict_7['codigo distrito'].append(l[13:14])
            dict_7['nombre territorio'].append(l[14:64].strip())
            dict_7['poblacion'].append(l[64:72])
            dict_7['numero mesas'].append(l[72:77])
            dict_7['censo ine'].append(l[77:85])
            dict_7['censo escrutinio'].append(l[85:93])
            dict_7['censo cere'].append(l[93:101])
            dict_7['votantes cere'].append(l[101:109])
            dict_7['votances primer avance'].append(l[109:117])
            dict_7['votantes segundo avance'].append(l[117:125])
            dict_7['votos blancos'].append(l[125:133])
            dict_7['votos nulos'].append(l[133:141])
            dict_7['votos candidaturas'].append(l[141:149])
            dict_7['numero escaños'].append(l[149:155])
            dict_7['votos afirmativos referendum'].append(l[155:163])
            dict_7['votos negativos referendum'].append(l[163:171])
            dict_7['datos oficiales'].append(l[171:172])

    df_7 = pd.DataFrame(dict_7)
    for col in range(len(dict_7_types)):
        col_name = df_7.columns[col]
        df_7[col_name] = df_7[col_name].astype(dict_7_types[col], errors='ignore')

    df_7.to_csv(carpeta + '07' + sufijo + '.csv', index=None, sep=';')

    # 8. CANDIDATURAS AMBITO SUPERIOR MUNICIPIO
    print('Procesando archivo 8')
    dict_8 = {'tipo eleccion': [], 'año eleccion': [], 'mes eleccion': [], 'numero vuelta': [],
              'codigo comunidad autonoma': [],
              'codigo provincia': [], 'numero distrito': [], 'codigo candidatura': [],
              'votos candidatura': [], 'numero candidatos': []}
    dict_8_types = [int] * 10

    with open(carpeta + '08' + sufijo + '.DAT', 'r', encoding="latin-1") as f:
        for l in f.readlines():
            dict_8['tipo eleccion'].append(l[:2])
            dict_8['año eleccion'].append(l[2:6])
            dict_8['mes eleccion'].append(l[6:8])
            dict_8['numero vuelta'].append(l[8:9])
            dict_8['codigo comunidad autonoma'].append(l[9:11])
            dict_8['codigo provincia'].append(l[11:13])
            dict_8['numero distrito'].append(l[13:14])
            dict_8['codigo candidatura'].append(l[14:20])
            dict_8['votos candidatura'].append(l[20:28])
            dict_8['numero candidatos'].append(l[28:33])

    df_8 = pd.DataFrame(dict_8)
    for col in range(len(dict_8_types)):
        col_name = df_8.columns[col]
        df_8[col_name] = df_8[col_name].astype(dict_8_types[col])

    df_8.to_csv(carpeta + '08' + sufijo + '.csv', index=None, sep=';')

    # 9.DATOS MESAS
    print('Procesando archivo 9')
    dict_9 = {'tipo eleccion': [], 'año eleccion': [], 'mes eleccion': [], 'numero vuelta': [],
              'codigo comunidad autonoma': [],
              'codigo provincia': [], 'codigo municipio': [], 'numero distrito': [], 'codigo seccion': [],
              'codigo mesa': [],
              'censo ine': [], 'censo escrutinio': [], 'censo cere': [], 'votantes cere': [],
              'votances primer avance': [], 'votantes segundo avance': [],
              'votos blancos': [], 'votos nulos': [], 'votos candidaturas': [], 'votos afirmativos referendum': [],
              'votos negativos referendum': [], 'datos oficiales': [], 'codigo': []}
    dict_9_types = [int] * 8 + [str] * 2 + [int] * 11 + 2 * [str]

    with open(carpeta + '09' + sufijo + '.DAT', 'r', encoding="latin-1") as f:
        for l in f.readlines():
            dict_9['tipo eleccion'].append(l[:2])
            dict_9['año eleccion'].append(l[2:6])
            dict_9['mes eleccion'].append(l[6:8])
            dict_9['numero vuelta'].append(l[8:9])
            dict_9['codigo comunidad autonoma'].append(l[9:11])
            dict_9['codigo provincia'].append(l[11:13])
            dict_9['codigo municipio'].append(l[13:16])
            dict_9['numero distrito'].append(l[16:18])
            dict_9['codigo seccion'].append(l[18:22])
            dict_9['codigo mesa'].append(l[22:23])
            dict_9['censo ine'].append(l[23:30])
            dict_9['censo escrutinio'].append(l[30:37])
            dict_9['censo cere'].append(l[37:44])
            dict_9['votantes cere'].append(l[44:51])
            dict_9['votances primer avance'].append(l[51:58])
            dict_9['votantes segundo avance'].append(l[58:65])
            dict_9['votos blancos'].append(l[65:72])
            dict_9['votos nulos'].append(l[72:79])
            dict_9['votos candidaturas'].append(l[79:86])
            dict_9['votos afirmativos referendum'].append(l[86:93])
            dict_9['votos negativos referendum'].append(l[93:100])
            dict_9['datos oficiales'].append(l[100:101])
            dict_9['codigo'].append(l[11:22].strip())

    df_9 = pd.DataFrame(dict_9)
    for col in range(len(dict_9_types)):
        col_name = df_9.columns[col]
        df_9[col_name] = df_9[col_name].astype(dict_9_types[col])

    df_9.to_csv(carpeta + '09' + sufijo + '.csv', index=None, sep=';')

    # 10 DATOS CANDIDATURAS MESA
    print('Procesando archivo 10')
    dict_10 = {'tipo eleccion': [], 'año eleccion': [], 'mes eleccion': [], 'numero vuelta': [],
               'codigo comunidad autonoma': [],
               'codigo provincia': [], 'codigo municipio': [], 'numero distrito': [], 'codigo seccion': [],
               'codigo mesa': [],
               'codigo candidatura': [], 'votos candidatura': [], 'codigo': []}
    dict_10_types = [int] * 8 + [str] * 2 + [int] * 2 + [str]

    with open(carpeta + '10' + sufijo + '.DAT', 'r', encoding="latin-1") as f:
        for l in tqdm(f.readlines()):
            dict_10['tipo eleccion'].append(l[:2])
            dict_10['año eleccion'].append(l[2:6])
            dict_10['mes eleccion'].append(l[6:8])
            dict_10['numero vuelta'].append(l[8:9])
            dict_10['codigo comunidad autonoma'].append(l[9:11])
            dict_10['codigo provincia'].append(l[11:13])
            dict_10['codigo municipio'].append(l[13:16])
            dict_10['numero distrito'].append(l[16:18])
            dict_10['codigo seccion'].append(l[18:22])
            dict_10['codigo mesa'].append(l[22:23])
            dict_10['codigo candidatura'].append(l[23:29])
            dict_10['votos candidatura'].append(l[29:36])
            dict_10['codigo'].append(l[11:22].strip())

    df_10 = pd.DataFrame(dict_10)
    for col in range(len(df_10.columns)):
        col_name = df_10.columns[col]
        df_10[col_name] = df_10[col_name].astype(dict_10_types[col])

    df_10.to_csv(carpeta + '10' + sufijo + '.csv', index=None, sep=';')

    """
    A partir de aquí vamos a hacer unas tablas más user-friendly. Para ello, vamos a crear una tabla con las secciones
    (que es lo que nos interesa), y vamos a incluir el recuento de votos tal cual lo obtenemos de la matriz de votos
    de candidatura. La matriz final tendrá el número de votantes (por CERE también), votos blancos, votos nulos, y resultados
    de los principales partidos políticos. La elección de los partidos la vamos a hacer a partir de un algoritmo genérico
    que describiré más abajo.
    """

    """
    Primero vamos a obtener los nombres de las candidaturas. Para ello vamos a crear un diccionario que transforme el número
    de la candidatura en las siglas. Debido a traducciones, existen partidos con el mismo nombre pero distinto número de candidatura.
    Por ejemplo, PODEMOS está en las formas PODEMOS-EUI, PODEMOS-EU, PODEMOS-IU, PODEMOS-IX y, además, PODEMOS-IU tiene 5
    candidaturas diferentes (3 unidas podemos, uno de alto Aragón y otra de Euskadi). Todo esto debería ser lo mismo. Así,
    deberíamos obtener una lista con las siglas iguales.

    Sin embargo, esto se complica aún más. PSOE por ejemplo tiene las formas PSOE (catalán y castellano), PSC, PSE, PSdeG-PSOE, etc.

    Por desgracia, no todos los partidos se pueden simplificar automáticamente, así que va a ir que ajustándolos manualmente.
    Con nuevas elecciones y nuevos partidos habrá que ir alterando las normas de igualdad de siglas.
    """

    df_candidaturas = df_3.copy()
    df_candidaturas = df_candidaturas[['codigo candidatura', 'siglas candidatura']]

    intercambios = {'AVANT ADELA': 'AVANT', 'AVANT-LOS V': 'AVANT', 'EB': 'EB', 'EB-AZ': 'EB',
                    'MÁS PAÍS-': 'MÁS PAÍS', 'M PAÍS': 'MÁS PAÍS', 'M PAÍS-': 'MÁS PAÍS',
                    'PCOE': 'PC', 'PCPA': 'PC', 'PCPC': 'PC', 'PCPE': 'PC', 'PCTE/ELAK': 'PC', 'PCTG': 'PC',
                    'UNIDOS PODEMOS': 'PODEMOS', 'UNIDAS PODEMOS': 'PODEMOS',
                    'PSC': 'PSOE', 'PSdeG-PSOE': 'PSOE', 'PSE-EE (PSO': 'PSOE', }

    lista_candts = [i.strip() for i in df_candidaturas['siglas candidatura'].values]
    lista_candts = [intercambios[i] if i in intercambios.keys() else i for i in lista_candts]

    # Reemplazamos todos los PODEMOS-XXXX y simplificamos siglas (las que se lleven escaño al menos). Incluimos, por facilidad, ECP-GUANYEM EL CAMBI por UP
    dict_siglas = {'PODEMOS-': 'UP', 'UNID': 'UP', 'ECP': 'UP', 'ERC-': 'ERC', 'JxCAT-': 'JxCAT', 
                   'CCa-': 'CCa', 'CC-': 'CCa', 'MÁS PAÍS': 'MP', 
                   'COMPROM': 'COMPR', 'CUP-': 'CUP', 'BNG-': 'BNG', '¡TERUEL': '¡TE!', 'PP-': 'PP'}

    for key, val in dict_siglas.items():
        lista_candts = [val if i.startswith(key) else i for i in lista_candts]

    df_candidaturas['siglas candidatura'] = lista_candts
    df_candidaturas['codigo candidatura'] = df_candidaturas['codigo candidatura'].astype(int)
    dict_siglas = dict(zip(df_candidaturas['codigo candidatura'].values, df_candidaturas['siglas candidatura'].values))

    """
    Ahora que tenemos las candidaturas, vamos a abrir el fichero de votos por candidatura y por sección.
    El objetivo final es crear una tabla por sección que recoja los votos por candidatura, los votos nulos y los blancos.
    El fichero de votos de candidatura no incluye los votos blancos y nulos. Estos los tendremos que extraer del fichero
    de los datos por mesa.
    Además, como por cada sección hay varias mesas, tenemos que sumar los votos de todas las mesas de esa sección.
    """

    # Primero sumamos los votos de candidaturas y mesas en general, para reducir el número de filas.
    df_votos_candidatura = df_10.copy()

    df_votos_candidatura = df_votos_candidatura[
        ['codigo', 'codigo candidatura', 'votos candidatura']]
    df_votos_candidatura = df_votos_candidatura.sort_values(
        by=['codigo', 'codigo candidatura', 'votos candidatura'])
    df_votos_candidatura = df_votos_candidatura[df_votos_candidatura['votos candidatura'] > 0]

    df_votos_candidatura = df_votos_candidatura.replace({'codigo candidatura': dict_siglas})
    df_votos_candidatura['votos candidatura'] = df_votos_candidatura.groupby(['codigo', 'codigo candidatura', ])[
        'votos candidatura'].transform('sum')
    df_votos_candidatura = df_votos_candidatura.drop_duplicates(keep='first')
    df_votos_candidatura = df_votos_candidatura.reset_index(drop=True)

    df_votos_candidatura = df_votos_candidatura.pivot(index='codigo', columns='codigo candidatura',
                                                      values='votos candidatura')

    df_votos_mesa = df_9.copy()
    df_votos_mesa = df_votos_mesa[
        ['codigo', 'censo escrutinio', 'votos blancos', 'votos nulos', 'votos candidaturas']]
    df_votos_mesa = df_votos_mesa.sort_values(
        by=['codigo', 'censo escrutinio', ])
    df_votos_mesa[['censo escrutinio', 'votos blancos', 'votos nulos', 'votos candidaturas']] = \
        df_votos_mesa.groupby(['codigo', ])[['censo escrutinio', 'votos blancos', 'votos nulos',
                                             'votos candidaturas']].transform('sum')

    df_votos_mesa = df_votos_mesa.drop_duplicates(keep='first')
    df_votos_mesa = df_votos_mesa.set_index('codigo')

    df_combinado = df_votos_candidatura.join(df_votos_mesa, how='outer', sort=True)

    df_combinado = df_combinado[df_votos_mesa.columns.tolist() + df_votos_candidatura.columns.tolist()]
    df_combinado = df_combinado.fillna(0).astype(int)
    
    # Este datframe está compuesto por todas las mesas electorales (~36k) y por los diferentes partidos (~70)
    # Para simplificar el guardado y acceso de datos, lo pasamos a sparse (si podemos) y guardamos como pickle.

    df_combinado.to_pickle(carpeta + '/resultados_candidaturas.pickle')


