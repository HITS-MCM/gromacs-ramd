from dash import Dash, html, dcc, callback, Output, Input
import plotly.express as px
import pandas as pd
import numpy as np

app = Dash(__name__)

app.layout = html.Div([
    html.H2(children='RAMD Dashboard', style={'textAlign':'center'}),
    dcc.Graph(id='graph-content'),
    dcc.Interval(
        id='interval-component',
        interval=5 * 1000, # in milliseconds
        n_intervals=0
    )
])

@callback(
    Output('graph-content', 'figure'),
    Input('interval-component', 'n_intervals')
)
def update_graph(value):
    ramd = np.loadtxt('build/tests/ramd-2.xvg', comments=['#', '@'])
    columns = ['time']
    for i in range(ramd[0].size - 1):
        columns.append('distance' + str(i))
    df = pd.DataFrame(ramd, columns=columns)
    return px.line(df, x='time', y=df.columns)

if __name__ == '__main__':
    app.run(debug=True, port=8051)
