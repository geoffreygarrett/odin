import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

ACCELERATION_FIELD = '/grid_points/accelerations'
POTENTIAL_FIELD = '/grid_points/potentials'
POSITION_FIELD = '/grid_points/positions'
ATTRIBUTE_FIELD = '/grid_points/positions'


def read_and_process_data(file_name):
    with h5py.File(file_name, 'r') as f:
        positions = np.array(f[POSITION_FIELD]).reshape(-1, 3)
        potentials = np.array(f[POTENTIAL_FIELD])
        accelerations = np.array(f[ACCELERATION_FIELD]).reshape(-1, 3)
        acceleration_magnitudes = np.linalg.norm(accelerations, axis=1)

        df = pd.DataFrame({
            'x_coordinate': positions[:, 0],
            'y_coordinate': positions[:, 1],
            'potential': potentials,
            'acceleration': list(accelerations),
            'acceleration_magnitude': acceleration_magnitudes
        })

        attrs = dict(f[ATTRIBUTE_FIELD].attrs)
        return df, attrs


def plot_field(df, field, ax, title):
    data = df.pivot(index='y_coordinate', columns='x_coordinate', values=field)
    contour = ax.contourf(data.columns, data.index, data.values, cmap='viridis')
    ax.set_title(title, fontsize=14)
    ax.set_xlabel('$x$ (m)', fontsize=12)
    ax.set_ylabel('$y$ (m)', fontsize=12)
    plt.colorbar(contour, ax=ax, format='%.2e')
    return ax


def plot_ellipse(ax, a, b):
    ellipse = Ellipse((0, 0), 2 * a, 2 * b, edgecolor='black', facecolor='orange', alpha=0.3,
                      linewidth=2, linestyle='dashed')
    ax.add_patch(ellipse)
    ax.set_aspect('equal')
    return ax


import sys

if __name__ == "__main__":
    # file_name = sys.argv[1]
    file_name = 'point_mass.h5'
    df, attrs = read_and_process_data(file_name)

    # Sort dataframe by x and y coordinates
    df = df.sort_values(['x_coordinate', 'y_coordinate'])

    # Ellipsoid parameters for Kleopatra
    a = attrs['a']  # Semi-major axis
    b = attrs['b']  # Semi-minor axis

    fig, axs = plt.subplots(1, 2, figsize=(15, 7))



    axs[0] = plot_field(df, 'potential', axs[0], 'Potential Field (XY Plane)')
    axs[0] = plot_ellipse(axs[0], a, b)

    axs[1] = plot_field(df, 'acceleration_magnitude', axs[1], 'Acceleration Field (XY Plane)')
    axs[1] = plot_ellipse(axs[1], a, b)

    # scatter plot grid points
    axs[0].scatter(df['x_coordinate'], df['y_coordinate'], c='black', s=1)

    plt.tight_layout()
    plt.show()
