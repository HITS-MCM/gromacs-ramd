from dash import Dash, html, dcc, callback, Output, Input
import pandas as pd
import numpy as np
import plotly.subplots as ps

app = Dash(__name__)

app.layout = html.Div(
    [
        html.H2(children="RAMD Dashboard", style={"textAlign": "center"}),
        dcc.Graph(id="graph-content", style={"width": "90vw", "height": "90vh"}),
        dcc.Interval(
            id="interval-component", interval=5 * 1000, n_intervals=0  # in milliseconds
        ),
    ]
)


@callback(Output("graph-content", "figure"), Input("interval-component", "n_intervals"))
def update_graph(value):
    ramd = np.loadtxt("build/tests/ramd-2.xvg", comments=["#", "@"])
    columns = ["time"]
    nb_dist = ramd[0].size - 1
    for i in range(nb_dist):
        columns.append("distance" + str(i))
    df = pd.DataFrame(ramd, columns=columns)

    fig = ps.make_subplots(
        rows=nb_dist,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.009,
        horizontal_spacing=0.009,
    )
    fig["layout"]["margin"] = {"l": 30, "r": 10, "b": 50, "t": 25}

    for r in range(nb_dist):
        fig.append_trace(
            {"x": df.time, "y": df.iloc[:, r + 1], "type": "scatter"}, r + 1, 1
        )

    return fig


if __name__ == "__main__":
    app.run(debug=True, port=8051)
