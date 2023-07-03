import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import pandas as pd
import numpy as np
import h5py

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


def read_trajectory_data(file_name, trajectory_number):
    with h5py.File(file_name, 'r') as f:
        trajectory_group = f[f'/trajectories/trajectory_{trajectory_number}']

        positions = np.array(trajectory_group['positions'])
        velocities = np.array(trajectory_group['velocities'])
        times = np.array(trajectory_group['times'])

        df = pd.DataFrame({
            'x_position': positions[:, 0],
            'y_position': positions[:, 1],
            'z_position': positions[:, 2],
            'x_velocity': velocities[:, 0],
            'y_velocity': velocities[:, 1],
            'z_velocity': velocities[:, 2],
            'time': times
        })

        attrs = dict(f['/trajectories'].attrs)
        return df, attrs


def plot_trajectory(df, ax):
    ax.plot(df['x_position'], df['y_position'], label='Spacecraft trajectory')
    ax.set_title('Spacecraft trajectory', fontsize=14)
    ax.set_xlabel('$x$ (m)', fontsize=12)
    ax.set_ylabel('$y$ (m)', fontsize=12)
    ax.legend()
    return ax


def plot_field(df, field, ax, title):
    vmax = df[field].max()
    vmin = df[field].min()
    data = df.pivot(index='y_coordinate', columns='x_coordinate', values=field)
    contour = ax.contourf(data.columns, data.index, data.values, cmap='viridis', vmin=vmin,
                          vmax=vmax)
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


if __name__ == "__main__":
    file_name1 = "point_mass.h5"
    file_name2 = "ellipsoid.h5"

    df1, attrs1 = read_and_process_data(file_name1)
    df2, attrs2 = read_and_process_data(file_name2)

    # Compute the differences
    df_diff = df1.copy()
    df_diff['potential'] = df1['potential'] - df2['potential']
    df_diff['acceleration_magnitude'] = df1['acceleration_magnitude'] - df2[
        'acceleration_magnitude']

    # Sort dataframes by x and y coordinates
    df1 = df1.sort_values(['x_coordinate', 'y_coordinate'])
    df2 = df2.sort_values(['x_coordinate', 'y_coordinate'])
    df_diff = df_diff.sort_values(['x_coordinate', 'y_coordinate'])

    # Ellipsoid parameters for Kleopatra
    a1, b1 = attrs1['a'], attrs1['b']  # Semi-major and minor axis for the first file
    a2, b2 = attrs2['a'], attrs2['b']  # Semi-major and minor axis for the second file

    fig, axs = plt.subplots(3, 2, figsize=(15, 21))

    # Plot the fields for the first file
    plot_field(df1, 'potential', axs[0, 0], f'Potential Field (XY Plane) - {file_name1}')
    plot_ellipse(axs[0, 0], a1, b1)
    plot_field(df1, 'acceleration_magnitude', axs[0, 1],
               f'Acceleration Field (XY Plane) - {file_name1}')
    plot_ellipse(axs[0, 1], a1, b1)

    # Plot the fields for the second file
    plot_field(df2, 'potential', axs[1, 0], f'Potential Field (XY Plane) - {file_name2}')
    plot_ellipse(axs[1, 0], a2, b2)
    plot_field(df2, 'acceleration_magnitude', axs[1, 1],
               f'Acceleration Field (XY Plane) - {file_name2}')
    plot_ellipse(axs[1, 1], a2, b2)

    # Plot the differences
    plot_field(df_diff, 'potential', axs[2, 0], 'Difference in Potential Field (XY Plane)')
    plot_field(df_diff, 'acceleration_magnitude', axs[2, 1],
               'Difference in Acceleration Field (XY Plane)')

    plt.tight_layout()
    plt.show()
