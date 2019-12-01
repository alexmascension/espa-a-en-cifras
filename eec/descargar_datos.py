"""
Este documento nos va a servir para descargar todos los datos basados en webs de descargas oficiales.
"""
import os
import urllib.request
import zipfile


def sacar_datos(directorio: str, fichero: str, web: str):
    if not os.path.exists(directorio):
        os.mkdir(directorio)

    urllib.request.urlretrieve(web, directorio + '/' + fichero)

    if fichero.endswith('.zip'):
        with zipfile.ZipFile(directorio + '/' + fichero, 'r') as z:
            z.extractall(directorio)
