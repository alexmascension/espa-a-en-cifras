import panel
import param
from bokeh.io import curdoc
from bokeh.themes import Theme

curdoc().theme = Theme(json={'attrs': {# apply defaults to Figure properties
'Figure': {
    'toolbar_location': 'right',
    'outline_line_color': None,
    'min_border_right': 10,
    'sizing_mode': 'stretch_width',
    'outline_line_color': None
},'Grid': {
    'grid_line_color': None,
},
'Title': {
    'text_font_size': '14pt'
},# apply defaults to Axis properties
'Axis': {
    'visible': False,
},
# apply defaults to Legend properties
'Legend': {
    'background_fill_alpha': 0.8,
}}})

def scatter_multiselect_subordinada(df, var_indepe, var_subord, x, y, titulo='', hover=None, alpha_max=0.9,
                                    alpha_min=0.02, **scatter_kwargs):
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

        def crear_plot(self, sub_df):
            scatter = figure(plot_height=400, plot_width=400, tools='reset,box_zoom,pan,wheel_zoom,lasso_select,undo,redo',
                             sizing_mode='scale_width', output_backend="webgl")
            scatter.add_tools(hover)
            scatter.scatter(x=x, y=y, source=sub_df, color='color', alpha='alpha', **scatter_kwargs)
            return scatter

        @panel.depends('indepe_obj', watch=True)
        def actualizar_subordinada(self):
            sub_df = df[df[var_indepe].isin(self.indepe_obj)]
            lista_subordinada = list(dict.fromkeys(sub_df[var_subord].values))
            self.param['subordinada_obj'].objects = lista_subordinada
            self.param['subordinada_obj'].default = lista_subordinada

        @panel.depends('indepe_obj', 'subordinada_obj')
        def plotear_auton(self):
            sub_df = df[df[var_indepe].isin(self.indepe_obj)]

            sub_df['alpha'] = alpha_min
            sub_df.loc[sub_df[var_subord].isin(self.subordinada_obj), 'alpha'] = alpha_max

            return self.crear_plot(sub_df)

    sss = ScatterAutonPart(name=titulo)
    panel_row = panel.Row(panel.Param(sss.param, widgets={'indepe_obj': {'size': 19, 'name': var_indepe},
                                     'subordinada_obj': {'size': 10, 'name': var_subord}
                                     }), sss.plotear_auton)
    return panel_row