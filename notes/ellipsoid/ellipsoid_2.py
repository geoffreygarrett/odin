import numpy as np
import plotly.graph_objs as go
import dash
from dash import dcc, html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc

def C(x, y, z, kappa, a, b, c):
    """Function to calculate C"""
    return (
            (x * x) / (a * a + kappa)
            + (y * y) / (b * b + kappa)
            + (z * z) / (c * c + kappa)
            - 1
    )

def calculate_segments(x, y, z, a, b, c):
    """Function to calculate C values for kappa in different segments"""
    segments = [
        np.linspace(-(c**2), np.max([x * a, y * b, z * c]), 1000),
        np.linspace(-(c**2), -(b**2), 1000),
        np.linspace(-(b**2), -(a**2), 1000),
        np.linspace(-(a**2), -np.max([x * a, y * b, z * c]) ** 2.5, 1000),
    ]
    C_segment_values = []
    for segment in segments:
        C_values = C(x, y, z, segment, a, b, c)
        C_segment_values.append(C_values)
    return segments, C_segment_values

app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

app.layout = dbc.Container([
    dcc.Graph(id='live-graph', animate=True),
    html.H2('Interactive Visualization of C(kappa) vs kappa', style={'textAlign': 'center'}),
    dbc.Row([
        dbc.Col([
            html.H4('Ellipsoid Parameters'),
            html.Label('a:'),
            dcc.Slider(
                id='slider-a',
                min=100,
                max=500,
                value=300,
                step=1,
                marks={i: str(i) for i in range(100, 501, 100)},
            ),
            html.Div(id='slider-a-display'),
            html.Label('b:'),
            dcc.Slider(
                id='slider-b',
                min=100,
                max=500,
                value=200,
                step=1,
                marks={i: str(i) for i in range(100, 501, 100)},
            ),
            html.Div(id='slider-b-display'),
            html.Label('c:'),
            dcc.Slider(
                id='slider-c',
                min=100,
                max=500,
                value=100,
                step=1,
                marks={i: str(i) for i in range(100, 501, 100)},
            ),
            html.Div(id='slider-c-display'),
            html.H4('Point Coordinates'),
            html.Label('x:'),
            dcc.Slider(
                id='slider-x',
                min=-500,
                max=500,
                value=-10,
                step=1,
                marks={i: str(i) for i in range(-500, 501, 100)},
            ),
            html.Div(id='slider-x-display'),
            html.Label('y:'),
            dcc.Slider(
                id='slider-y',
                min=-500,
                max=500,
                value=-10,
                step=1,
                marks={i: str(i) for i in range(-500, 501, 100)},
            ),
            html.Div(id='slider-y-display'),
            html.Label('z:'),
            dcc.Slider(
                id='slider-z',
                min=-500,
                max=500,
                value=-10,
                step=1,
                marks={i: str(i) for i in range(-500, 501, 100)},
            ),
            html.Div(id='slider-z-display'),
        ], md=4),
        dbc.Col([
            dcc.Graph(id='live-graph', config={'displayModeBar': False})
        ], md=8)
    ])
], fluid=True)

@app.callback(
    [Output('live-graph', 'figure')] +
    [Output(f'slider-{i}-display', 'children') for i in 'abcxyz'],
    [Input(f'slider-{i}', 'value') for i in 'abcxyz']
)
def update_graph(a, b, c, x, y, z):
    kappa_segments, C_value_segments = calculate_segments(x, y, z, a, b, c)

    fig = go.Figure()

    colors = ['blue', 'green', 'red', 'purple']
    labels = ['Ellipsoidal Confocal Family', 'Hyperboloid 1', 'Hyperboloid 2', 'Fourth segment']

    for kappa_values, C_values, color, label in zip(kappa_segments, C_value_segments, colors, labels):
        fig.add_trace(go.Scatter(x=kappa_values, y=C_values, mode='lines', name=label, line=dict(color=color)))
        crossing_indices = np.where(np.diff(np.sign(C_values)))[0]
        roots = kappa_values[crossing_indices]
        fig.add_trace(go.Scatter(x=roots, y=[0]*len(roots), mode='markers', name='Roots', marker=dict(color='red')))

    fig.update_layout(height=600, width=800, title_text="Plot of C(kappa) vs kappa")

    return [fig] + [f"Current value: {val}" for val in [a, b, c, x, y, z]]

if __name__ == '__main__':
    app.run_server(debug=True)
