import panel
import param
from bokeh.models import HoverTool, LinearColorMapper, LogColorMapper
from bokeh.io import show, output_file, push_notebook, curdoc
from bokeh.models import ColumnDataSource, Select, HoverTool, Panel, Tabs, Label, ColorBar, BasicTicker, PrintfTickFormatter
from bokeh.plotting import figure
from bokeh.models.widgets import MultiSelect, CheckboxGroup, Dropdown
from bokeh.layouts import column, row
from bokeh.themes import Theme
from bokeh.models import Range1d

import pandas as pd
import numpy as np

from scipy.stats import linregress as lr
from scipy.optimize import curve_fit as cf


def scatter_multiselect_subordinada(df, var_indepe, var_subord, x, y, titulo='', alpha_max=0.9, corr_type='lineal',
                                    alpha_min=0.04, **scatter_kwargs):
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
            y_start, y_end = min(sub_df[y][sub_df['alpha'] == alpha_max]), max(sub_df[y][sub_df['alpha'] == alpha_max])
            y_start -= 0.05 * (y_end - y_start)
            y_end += 0.05 * (y_end - y_start)

            scatter = figure(plot_height=400, plot_width=400,
                             tools='reset,box_zoom,pan,wheel_zoom,lasso_select,undo,redo',
                             sizing_mode='scale_width', output_backend="webgl",
                             y_range=(y_start, y_end))

            hover = HoverTool(
                tooltips="""
                    <div><span style="font-size: 17px; font-weight: bold;">@municipio</span></div>
                    <div><span style="font-size: 12px;">@provincia (@autonomia), @poblacion</span></div>
                    <div><span style="font-size: 14px; font-weight: bold;">@partido</span>
                    <span style="font-size: 12px;">@porcentaje</span></div>
                    """)

            scatter.add_tools(hover)
            scatter.scatter(x=x, y=y, source=sub_df, color='color', alpha='alpha', **scatter_kwargs)
            y_range=(y_start, y_end)
            if corr:  # Con esto añadimos las correlaciones.
                for var_subor_i in self.subordinada_obj:
                    for var_indep_i in self.indepe_obj:
                        si_df = sub_df[(sub_df[var_indepe] == var_indep_i) &
                                       (sub_df[var_subord] == var_subor_i)]

                        if len(si_df) > 1:  # Nos aseguramos porque si no falla para Ceuta y Melilla
                            x_vals, y_vals = si_df[x].values, si_df[y].values

                            if corr_type == 'lineal':
                                def f(x, m, b):
                                    return m * x + b

                                m, b, r, p, err = lr(x_vals, y_vals)
                                text_label = "r² = %.2f" % r**2
                            elif corr_type == 'exp':
                                def f(x, m, b):
                                    return np.power(m * np.log10(x) + b, 10)

                                popt, pcor = cf(f, x_vals, y_vals)
                                m, b = popt
                                ss_res = np.sum((y_vals - f(x_vals, m, b)) ** 2)
                                ss_tot = np.sum((y_vals - np.mean(y_vals)) ** 2)
                                r_squared = 1 - (ss_res / ss_tot)
                                text_label = "r² = %.2f" % r_squared

                            x_arr = np.linspace(min(x_vals), max(x_vals), 100)
                            scatter.line(x_arr, [f(x_i, m, b) for x_i in x_arr], color=si_df['color'].iloc[0])
                            r_label = Label(x=1.05 * max(x_vals), y=f(max(x_vals), m, b),
                                            text=text_label, text_align='left', text_color=si_df['color'].iloc[0],
                                            render_mode='css')
                            scatter.add_layout(r_label)
                            if f(max(x_vals), m, b) > y_end: y_end = f(max(x_vals), m, b)
                            if f(max(x_vals), m, b) < y_start: y_start = f(max(x_vals), m, b)
                        else:
                            pass
                scatter.y_range = Range1d(y_start, y_end)
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


import panel
import param
from bokeh.models import HoverTool, LinearColorMapper, LogColorMapper
from bokeh.io import show, output_file, push_notebook, curdoc
from bokeh.models import ColumnDataSource, Select, HoverTool, Panel, Tabs, Label, ColorBar, BasicTicker, \
    PrintfTickFormatter
from bokeh.plotting import figure
from bokeh.models.widgets import MultiSelect, CheckboxGroup, Dropdown
from bokeh.layouts import column, row
from bokeh.themes import Theme

import pandas as pd
import numpy as np

from scipy.stats import linregress as lr
from scipy.optimize import curve_fit as cf


def heatmap_correlaciones(df, var_x_ht, var_y_ht, var_x, var_y, cmap='div', corr_type='lineal'):
    x, y, z = [], [], []
    opts_x, opts_y = sorted(list(dict.fromkeys(df[var_x_ht].values))), sorted(list(dict.fromkeys(df[var_y_ht].values)))[
                                                                       ::-1]
    #     opts_x = ['PSOE', 'PP']
    #     opts_y = ['Euskadi', 'Cataluña']
    for opt_x in opts_x:
        for opt_y in opts_y:
            df_xy = df[(df[var_x_ht] == opt_x) & (df[var_y_ht] == opt_y)]
            x.append(opt_x);
            y.append(opt_y)

            try:
                x_vals, y_vals = df_xy[var_x].values, df_xy[var_y].values

                if corr_type == 'lineal':
                    m, b, r, p, err = lr(x_vals, y_vals)
                    z.append(r ** 2)
                elif corr_type == 'exp':
                    def f(x, m, b):
                        return np.power(m * np.log10(x) + b, 10)

                    popt, pcor = cf(f, x_vals, y_vals)
                    m, b = popt
                    ss_res = np.sum((y_vals - f(x_vals, m, b)) ** 2)
                    ss_tot = np.sum((y_vals - np.mean(y_vals)) ** 2)
                    r_squared = 1 - (ss_res / ss_tot)
                    z.append(r_squared)
            except:
                z.append(np.NaN)

    df_correlacion = pd.DataFrame({var_x_ht: x, var_y_ht: y, 'z': z})
    p = figure(plot_height=150, plot_width=400, x_range=opts_x, y_range=opts_y,
               tools='hover',
               sizing_mode='scale_width', tooltips=[(var_x_ht, '@' + var_x_ht),
                                                    (var_y_ht, '@' + var_y_ht),
                                                    ('r', '@z')])

    if cmap == 'seq':
        colors = ['#d1eeea', '#a8dbd9', '#85c4c9', '#68abb8', '#4f90a6', '#3b738f', '#2a5674']
    elif cmap == 'div':
        colors = ['#009B9E', '#42B7B9', '#A7D3D4', '#F1F1F1', '#E4C1D9', '#D691C1', '#C75DAB']
    else:
        colors = cmap

    mapper = LinearColorMapper(palette=colors, low=df_correlacion.z.min(), high=df_correlacion.z.max(),
                               nan_color="#bbbbbb")

    p.grid.grid_line_color = None
    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "8pt"
    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = 0.8  # pi/4

    p.rect(source=df_correlacion, x=var_x_ht, y=var_y_ht, width=1, height=1,
           fill_color={'field': 'z', 'transform': mapper}, line_color=None)

    color_bar = ColorBar(color_mapper=mapper, major_label_text_font_size="5pt",
                         ticker=BasicTicker(desired_num_ticks=len(colors)),
                         formatter=PrintfTickFormatter(format="%.3f"),
                         label_standoff=6, border_line_color=None, location=(0, 0))
    p.add_layout(color_bar, 'right')

    return p