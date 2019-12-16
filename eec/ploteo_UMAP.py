import param
import panel
import holoviews as hv
import pandas as pd
import numpy as np
import panel
import param
from bokeh.io import curdoc
from bokeh.themes import Theme
from bokeh.plotting import figure
import matplotlib.pyplot as plt
import matplotlib as mpl

import umap
import scipy.stats as sts
import scanpy as sc
from bokeh.models import HoverTool

dict_prov = {'Araba/Álava': 'Euskadi', 'Albacete': "Castilla la Mancha", 'Alicante': "Comunidad Valenciana",'Almería': 'Andalucía',
 'Ávila': "Castilla y León", 'Badajoz': "Extremadura", 'Balears, Illes': "Baleares", 'Barcelona': "Cataluña",
 'Burgos': "Castilla y León", 'Cáceres': "Extremadura", 'Cádiz': 'Andalucía', 'Castellón/Castelló': "Comunidad Valenciana",
 'Ciudad Real': "Castilla la Mancha", 'Córdoba': 'Andalucía', 'Coruña, A': "Galicia",
 'Cuenca': "Castilla la Mancha", 'Girona': "Cataluña", 'Granada': 'Andalucía',
 'Guadalajara': "Castilla la Mancha", 'Gipuzkoa': 'Euskadi', 'Huelva': 'Andalucía', 'Huesca': 'Aragón',
 'Jaén': 'Andalucía', 'León': "Castilla y León", 'Lleida': "Cataluña",
 'Rioja, La': "La Rioja", 'Lugo': "Galicia", 'Madrid': "Madrid", 'Málaga': 'Andalucía',
 'Murcia': "Murcia", 'Navarra': "Navarra", 'Ourense': "Galicia", 'Asturias': "Asturias",
 'Palencia': "Castilla y León", 'Palmas, Las':'Canarias', 'Pontevedra': "Galicia", 'Salamanca': "Castilla y León",
 'Santa Cruz de Tenerife': 'Canarias', 'Cantabria': "Cantabria", 'Segovia': "Castilla y León", 'Sevilla': "Andalucía",
 'Soria': "Castilla y León", 'Tarragona': "Cataluña", 'Teruel': 'Aragón', 'Toledo': "Castilla la Mancha",
 'Valencia/Valéncia': "Comunidad Valenciana", 'Valladolid': "Castilla y León", 'Bizkaia': 'Euskadi',
 'Zamora': "Castilla y León", 'Zaragoza': 'Aragón', 'Ceuta': 'Ceuta', 'Melilla': 'Melilla'}

dict_num_prov = {list(dict_prov.keys())[i]:i for i in range(len(dict_prov))}



def calcular_leiden(array, res, subres, seed):
    print('Calculando leiden')

    adata = sc.AnnData(X=np.nan_to_num(array))
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=res, random_state=seed)

    if subres > 0:
        array_return = np.zeros(len(adata))
        clusters = list(dict.fromkeys(adata.obs['leiden'].values))
        n_clusters = 0

        for cluster in clusters:
            index_cluster = np.argwhere(adata.obs['leiden'] == cluster).flatten()
            subadata = adata[index_cluster].copy()
            subadata.X = np.nan_to_num(subadata.X)
            sc.pp.neighbors(subadata)
            sc.tl.leiden(subadata, resolution=subres, random_state=seed)
            array_return[index_cluster] = subadata.obs['leiden'].values.astype(int) + n_clusters
            n_clusters += len(list(dict.fromkeys(subadata.obs['leiden'])))

        return array_return
    else:
        return adata.obs['leiden'].values.astype(int)


def calcular_UMAP(df, seed):
    print('Calculando UMAP')
    array_porcentajes = df.loc[:,
                        [i for i in df.columns if (i.startswith('porcentaje_')) & ((not i[-1].isdigit()))]].values
    array_porcentajes = np.nan_to_num(array_porcentajes)

    reducer = umap.UMAP(transform_seed=seed)
    embedding = reducer.fit_transform(array_porcentajes)

    return embedding


def map_sizes(r, min_pop, max_pop, r_min=2, r_max=13,):
    return r_min + (r - min_pop) / (max_pop - min_pop) * (r_max - r_min)


def crear_df_datos(df, resolucion, subresolucion, seed, r_min=2, r_max=13, ):
    embedding = calcular_UMAP(df, seed)
    array = df[[i for i in df.columns
                if (i.startswith('porcentaje_')) ]].copy().values
    df['leiden'] = calcular_leiden(array, resolucion, subresolucion, seed)
    min_pop, max_pop = np.power(min(df['censo escrutinio'].values), 0.05), np.power(max(df['censo escrutinio'].values),
                                                                                    0.05)
    df_datos = pd.DataFrame(dict({'x': embedding[:, 0], 'y': embedding[:, 1],
                                  'partido_1': df['partido_1'].values,
                                  'partido_2': df['partido_2'].values,
                                  'partido_3': df['partido_3'].values,
                                  'autonomia': df['Autonomia'], 'provincia': df['CPRO'],
                                  'color': ['#000000'] * len(df), 'Nombre': df['Nombre'],
                                  'poblacion': df['censo escrutinio'].values,
                                  'tamano': [map_sizes(np.power(i, 0.05), min_pop, max_pop, r_min, r_max)
                                             for i in df['censo escrutinio'].values],
                                  'alpha': [0.99 for _ in df['censo escrutinio'].values],
                                  'leiden': df['leiden'].values},
                                 **{i: df[i].values for i in df.columns if (i.startswith('porcentaje_'))}))

    return df_datos


def dibujar_UMAP_votos_autonomia(df, dict_colores, alpha_max=0.95, alpha_min=0.00, titulo=''):
    class UMAPPlotClass(param.Parameterized):
        lista_opciones = ['leiden'] + ['partido_1', 'partido_2', 'partido_3'] + ['provincia', 'autonomia'] + \
                         [i for i in df.columns if (i.startswith('porcentaje_')) &
                          (i not in ['porcentaje_1', 'porcentaje_2', 'porcentaje_3'])]
        lista_autonomias = sorted(list(dict.fromkeys(df['autonomia'].values)))

        obj_opciones = param.ObjectSelector(default='leiden', objects=lista_opciones)
        obj_autonomias = param.ListSelector(default=lista_autonomias, objects=lista_autonomias)

        def color_mapper(self, c2, c1, mix):
            if np.isnan(mix): mix = 0
            c1 = np.array(mpl.colors.to_rgb(c1))
            c2 = np.array(mpl.colors.to_rgb(c2))
            return mpl.colors.to_hex((1 - mix) * c1 + mix * c2)

        def crear_plot(self, sub_df):
            UMAP = figure(plot_height=700, plot_width=700, tools='box_zoom,reset,pan,wheel_zoom,lasso_select,undo,redo',
                          sizing_mode='scale_width', output_backend="webgl", toolbar_location='right')

            hover_UMAP = HoverTool(
                tooltips="""
                    <div><span style="font-size: 17px; font-weight: bold;">@Nombre</span></div>
                    <div><span style="font-size: 12px;">@provincia (@autonomia), @poblacion</span></div>
                    <div><span style="font-size: 14px; font-weight: bold;">@partido_1</span>
                    <span style="font-size: 13px;">@porcentaje_1 %</span></div>
                    <div><span style="font-size: 14px; font-weight: bold;">@partido_2</span>
                    <span style="font-size: 13px;">@porcentaje_2 %</span></div>
                    <div><span style="font-size: 14px; font-weight: bold;">@partido_3</span>
                    <span style="font-size: 13px;">@porcentaje_3 %</span></div>
                    <div><span style="font-size: 14px; font-weight: bold;">% Abstención</span>
                    <span style="font-size: 13px;">@porcentaje_abstencion %</span></div>
                    <div><span style="font-size: 14px; font-weight: bold;">Grupo</span>
                    <span style="font-size: 13px;">@leiden</span></div>""",
            )

            UMAP.add_tools(hover_UMAP)
            UMAP.axis.visible, UMAP.xgrid.visible, UMAP.ygrid.visible = False, False, False
            UMAP.scatter('x', 'y', source=sub_df, line_alpha='alpha_line', line_width=0.3, line_color="#000000",
                         size='tamano', color='color', alpha='alpha')

            return UMAP

        @panel.depends('obj_opciones', 'obj_autonomias')
        def update_plot(self):
            dict_color_leiden = {0: '#5F4690', 1: '#1D6996', 2: '#38A6A5', 3: '#0F8554', 4: '#73AF48',
                                 5: '#EDAD08', 6: '#E17C05', 7: '#CC503E', 8: '#94346E', 9: '#6F4070',
                                 10: '#994E95', 11: '#666666', 12: '#11A579', 13: '#E73F74', 14: '#CF1C90',
                                 15: '#3969AC', 16: '#4b4b8f', 17: '#66C5CC', 18: '#F89C74', 19: '#DCB0F2',
                                 20: '#FE88B1', 21: '#661100'}
            dict_num_autonomias = {'Navarra': 0, 'Canarias': 1, 'Baleares': 2, 'Comunidad Valenciana': 3,
                                   'Asturias': 4, 'Andalucía': 5, 'Aragón': 6, 'Murcia': 7, 'Extremadura': 8,
                                   'Ceuta': 9, 'Cataluña': 10,
                                   'Cantabria': 11, 'Madrid': 12, 'Castilla y León': 13, 'La Rioja': 14, 'Euskadi': 15,
                                   'Melilla': 16,
                                   'Castilla la Mancha': 17, 'Galicia': 18}
            sub_df = df.copy()
            # Primero aplicamos el color
            attr = self.obj_opciones
            if attr in ['partido_1', 'partido_2', 'partido_3']:
                sub_df['color'] = [dict_colores[i] for i in sub_df[attr]]
            elif attr in ['porcentaje_blancos', 'porcentaje_nulos', 'porcentaje_abstencion']:
                max_attr = max(sub_df[attr])
                sub_df['color'] = [self.color_mapper('#000000', "#bbbbbb", i / max_attr) for i in sub_df[attr]]
            elif attr == 'leiden':
                sub_df['color'] = [dict_color_leiden[i % len(dict_color_leiden)] for i in sub_df[attr]]
            elif attr == 'provincia':
                sub_df['color'] = [dict_color_leiden[dict_num_prov[i] % 15] for i in sub_df['provincia']]
            elif attr == 'autonomia':
                sub_df['color'] = [dict_color_leiden[dict_num_autonomias[i] % 19] for i in sub_df[attr]]
            else:
                max_attr = max(sub_df[attr])
                sub_df['color'] = [self.color_mapper(dict_colores[attr.replace('porcentaje_', '')],
                                                     "#f0f0f0", i / max_attr) for i in sub_df[attr]]
            # Ahora aplicamos el alpha de las autonomías
            aut_select = self.obj_autonomias
            sub_df['alpha'] = alpha_min
            sub_df['alpha_line'] = alpha_min
            sub_df.loc[sub_df['autonomia'].isin(aut_select), 'alpha'] = alpha_max
            sub_df.loc[sub_df['autonomia'].isin(aut_select), 'alpha_line'] = 0.45

            return self.crear_plot(sub_df)

    sss = UMAPPlotClass(name=titulo)
    panel_row = panel.Row(panel.Param(sss.param, widgets={'obj_opciones': {'name': 'Opciones'},
                                                          'obj_autonomias': {'size': 19, 'name': 'Autonomías'}
                                                          }), sss.update_plot)
    return panel_row