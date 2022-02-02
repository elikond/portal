from __future__ import generator_stop
from typing import Optional

import dash
from dash.dependencies import Input, Output
from dash.exceptions import PreventUpdate
import dash_core_components as dcc
import dash_html_components as html

import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import os
from collections import Counter
from plotly.subplots import make_subplots
from pathlib import Path

import heatmap
from heatmap import *

fig = go.Figure()

def configure_app(app: dash.Dash):
    class Ids:
        pass

    app.layout = html.Div(
        children = [
            html.Div(children=[
                html.Label('Virus 1'),
                dcc.Dropdown(
                    id="selected-vir",
                    options=[
                        {'label': 'Dengue', 'value': 'DENV'},
                        {'label': 'Enterovirus', 'value': 'EV'},
                        {'label': 'Hepatitis A', 'value': 'HAV'},
                        {'label': 'Hepatitis C', 'value': 'HCV'},
                        {'label': 'Rhinovirus', 'value': 'RV'},
                        {'label': 'Human Coronavirus 229E', 'value': 'Wang_229E'},
                        {'label': 'Human Coronavirus OC43', 'value': 'Wang_OC43'},
                        {'label': 'SARS-CoV-2', 'value': 'Wang_SARS-CoV2'},
                    ],
                ),
            ],
        ),
                    html.Div(children=[
                html.Label('Virus 2'),
                dcc.Dropdown(
                    id="selected-vir1",
                    options=[
                        {'label': 'Dengue', 'value': 'DENV'},
                        {'label': 'Enterovirus', 'value': 'EV'},
                        {'label': 'Hepatitis A', 'value': 'HAV'},
                        {'label': 'Hepatitis C', 'value': 'HCV'},
                        {'label': 'Rhinovirus', 'value': 'RV'},
                        {'label': 'Human Coronavirus 229E', 'value': 'Wang_229E'},
                        {'label': 'Human Coronavirus OC43', 'value': 'Wang_OC43'},
                        {'label': 'SARS-CoV-2', 'value': 'Wang_SARS-CoV2'},
                    ],
                ),
            ],
        ),
            html.Div(
                className="graphContainer",
                children=[
                    dcc.Graph(className="graph", id="heatmap"),
            ],
            )
        ],
        )

    @app.callback(
        Output("heatmap","figure"),
        Input("selected-vir","value"),
        Input("selected-vir1","value")
    )

    def update_figure(vir1, vir2):
        final_df, a, combined_df = final([vir1, vir2], tot_vir)
        if not final_df.empty:
                fig = px.imshow(final_df, labels=dict(x="Viruses", y="Genes", color="Significance (-log[pos score])"),
                y=combined_df['Shared_Genes'][a], x = [vir1, vir2], title=combined_df['Original Name_x'][a])
                return fig
        else:
            data = [go.Heatmap( x=[], y=[], z=[])]
            fig = go.Figure(data=data)

            fig.update_layout(
                title = 'No data to display'
            )
            return fig

    return app

def run():
    app = dash.Dash(__name__)

    configure_app(app)

    app.run_server(debug=True)