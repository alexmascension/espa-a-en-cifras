import panel
import param
from bokeh.models import HoverTool, LinearColorMapper, LogColorMapper
from bokeh.io import show, output_file, push_notebook, curdoc
from bokeh.models import ColumnDataSource, Select, HoverTool, Panel, Tabs, Label
from bokeh.plotting import figure
from bokeh.models.widgets import MultiSelect, CheckboxGroup, Dropdown
from bokeh.layouts import column , row
from bokeh.themes import Theme

import pandas as pd
import numpy as np

def scatter_multiselect_subordinada(df, var_indepe, var_subord, x, y, titulo='', alpha_max=0.9,
                                    alpha_min=0.05, **scatter_kwargs):
    from scipy.stats import linregress as lr
    """
    Esta función está diseñada para dibujar un scatter con dos variables. La primera es una variable independiente,
    que actualiza el scatter cada vez que se  actualiza. Cuando la primera variable se actializa, da a elegir varias
    opciones en la segunda, de ahí que se denomine "subordinada".
    Una vez la variable subordinada ha sido actualizada, ésta se puede cambiar, que a su vez altera el scatter.
    En este caso, la acción de la variable subordinada va a ser alterar el valor alpha del scatter. Entonces,
    vamos a establecer un alpha_max para los puntos en general, y un alpha_min para los puntos que no han sido
    seleccionados.
    El color tiene que estar almacenado en la columna "color" del dataframe.
    """

    class ScatterAutonPart(param.Parameterized):
        lista_indepe = sorted(list(dict.fromkeys(df[var_indepe].values)))
        lista_subordinada = list(dict.fromkeys(df[var_subord].values))

        lista_opciones_subordinada = list(dict.fromkeys(df[var_subord].values))
        indepe_obj = param.ListSelector(default=lista_indepe, objects=lista_indepe)
        subordinada_obj = param.ListSelector(default=lista_subordinada, objects=lista_subordinada)
        correlacion_obj = param.Boolean(False, doc="Mostrar correlación")

        def crear_plot(self, sub_df, corr):
            scatter = figure(plot_height=400, plot_width=400,
                             tools='reset,box_zoom,pan,wheel_zoom,lasso_select,undo,redo',
                             sizing_mode='scale_width', output_backend="webgl")

            hover = HoverTool(
                tooltips="""
                    <div><span style="font-size: 17px; font-weight: bold;">@municipio</span></div>
                    <div><span style="font-size: 12px;">@provincia (@autonomia), @poblacion</span></div>
                    <div><span style="font-size: 14px; font-weight: bold;">@partido</span>
                    <span style="font-size: 12px;">@porcentaje</span></div>
                    """)

            scatter.add_tools(hover)
            scatter.scatter(x=x, y=y, source=sub_df, color='color', alpha='alpha', **scatter_kwargs)

            if corr: # Con esto añadimos las correlaciones.
                for var_subor_i in self.subordinada_obj:
                    for var_indep_i in self.indepe_obj:
                        si_df = sub_df[(sub_df[var_indepe] == var_indep_i) &
                                       (sub_df[var_subord] == var_subor_i)]

                        if len(si_df) > 1: # Nos aseguramos porque si no falla para Ceuta y Melilla
                            x_vals, y_vals = si_df[x].values, si_df[y].values
                            m, b, r, p, err = lr(x_vals, y_vals)
                            scatter.line([min(x_vals), max(x_vals)], [m * min(x_vals) + b, m * max(x_vals) + b],
                                         color=si_df['color'].iloc[0])
                            r_label = Label(x=1.01 * max(x_vals), y=m * max(x_vals) + b,
                                            text="r² = ." + ('0' + str(int(100 * r ** 2)))[-2:],
                                            text_align='left', text_color=si_df['color'].iloc[0],
                                            angle=np.arctan(m * (max(sub_df[x].values) - min(sub_df[x].values)) / (
                                                    max(sub_df[y].values) - min(sub_df[y].values))),
                                            angle_units='rad', render_mode='css')
                            scatter.add_layout(r_label)
                        else:
                            pass
            return scatter

        @panel.depends('indepe_obj', watch=True)
        def actualizar_subordinada(self):
            sub_df = df[df[var_indepe].isin(self.indepe_obj)]
            lista_subordinada = list(dict.fromkeys(sub_df[var_subord].values))
            self.param['subordinada_obj'].objects = lista_subordinada
            self.param['subordinada_obj'].default = lista_subordinada

        @panel.depends('indepe_obj', 'subordinada_obj', 'correlacion_obj')
        def plotear_auton(self):
            sub_df = df[df[var_indepe].isin(self.indepe_obj)]

            sub_df['alpha'] = alpha_min
            sub_df.loc[sub_df[var_subord].isin(self.subordinada_obj), 'alpha'] = alpha_max

            si = self.correlacion_obj
            return self.crear_plot(sub_df, si)

    sss = ScatterAutonPart(name=titulo)
    panel_row = panel.Row(panel.Param(sss.param, widgets={
        'indepe_obj': {'size': 19, 'name': var_indepe},
        'subordinada_obj': {'size': 10, 'name': var_subord},
        'correlacion_obj': {'name': 'Mostrar correlación'}
    }), sss.plotear_auton)

    return panel_row