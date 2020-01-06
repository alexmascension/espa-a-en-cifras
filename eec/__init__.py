from .descargar_datos import sacar_datos
from .procesar_datos_demograficos import procesar_datos_2011
from .procesar_renta_demografia import procesar_demografia, procesar_renta
from .procesar_elecciones import procesar_elecciones
from .funciones_ploteo import scatter_multiselect_subordinada, heatmap_correlaciones
from .ploteo_UMAP import dibujar_UMAP_votos_autonomia, crear_df_datos
# from .procesar_datos_demograficos import pro

__author__ = ', '.join([
    'Alex M. Ascensión',
])
__email__ = ', '.join([
    'alexmascension@gmail.com',
    # We don’t need all, the main authors are sufficient.
])

__version__ = "1.0"
