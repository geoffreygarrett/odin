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


def plot_field(df1, df2, field, ax, title):
    data1 = df1.pivot(index='y_coordinate', columns='x_coordinate', values=field)
    data2 = df2.pivot(index='y_coordinate', columns='x_coordinate', values=field)
    contour1 = ax.contourf(data1.columns, data1.index, data1.values, cmap='viridis')
    contour2 = ax.contourf(data2.columns, data2.index, data2.values, cmap='viridis')
    ax.set_title(title, fontsize=14)
    ax.set_xlabel('$x$ (m)', fontsize=12)
    ax.set_ylabel('$y$ (m)', fontsize=12)
    plt.colorbar(contour1, ax=ax, format='%.2e')
    plt.colorbar(contour2, ax=ax, format='%.2e')
    return ax


def plot_ellipse(ax, a, b):
    ellipse = Ellipse((0, 0), 2 * a, 2 * b, edgecolor='black', facecolor='orange', alpha=0.3,
                      linewidth=2, linestyle='dashed')
    ax.add_patch(ellipse)
    ax.set_aspect('equal')
    return ax


import sys

if __name__ == "__main__":
    # get file names from command line arguments
    file_name2 = "point_mass.h5"
    file_name1 = "ellipsoid.h5"

    df1, attrs1 = read_and_process_data(file_name1)
    df2, attrs2 = read_and_process_data(file_name2)

    # Sort dataframes by x and y coordinates
    df1 = df1.sort_values(['x_coordinate', 'y_coordinate'])
    df2 = df2.sort_values(['x_coordinate', 'y_coordinate'])

    # Ellipsoid parameters for Kleopatra
    a1, b1 = attrs1['a'], attrs1['b']  # Semi-major and minor axis for the first file
    a2, b2 = attrs2['a'], attrs2['b']  # Semi-major and minor axis for the second file

    fig, axs = plt.subplots(2, 2, figsize=(15, 14))

    axs[0, 0] = plot_field(df1, df2, 'potential', axs[0, 0],
                           f'Potential Field (XY Plane) - {file_name1}')
    axs[0, 0] = plot_ellipse(axs[0, 0], a1, b1)
    axs[0, 1] = plot_field(df1, df2, 'acceleration_magnitude', axs[0, 1],
                           f'Acceleration Field (XY Plane) - {file_name1}')
    axs[0, 1] = plot_ellipse(axs[0, 1], a1, b1)
    axs[1, 0] = plot_field(df1, df2, 'potential', axs[1, 0],
                           f'Potential Field (XY Plane) - {file_name2}')
    axs[1, 0] = plot_ellipse(axs[1, 0], a2, b2)
    axs[1, 1] = plot_field(df1, df2, 'acceleration_magnitude', axs[1, 1],
                           f'Acceleration Field (XY Plane) - {file_name2}')
    axs[1, 1] = plot_ellipse(axs[1, 1], a2, b2)

    # # scatter plot grid points for both dataframes
    # axs[0, 0].scatter(df1['x_coordinate'], df1['y_coordinate'], c='black', s=1)
    # axs[0, 1].scatter(df1['x_coordinate'], df1['y_coordinate'], c='black', s=1)
    # axs[1, 0].scatter(df2['x_coordinate'], df2['y_coordinate'], c='black', s=1)
    # axs[1, 1].scatter(df2['x_coordinate'], df2['y_coordinate'], c='black', s=1)

    plt.tight_layout()
    plt.show()
