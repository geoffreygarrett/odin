import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


def read_and_process_data(file_name):
    df = pd.read_csv(file_name)
    df['acceleration_magnitude'] = np.sqrt(df['x_acceleration'] ** 2 + df['y_acceleration'] ** 2)
    return df


def plot_field(df, field, ax, title):
    data = df.pivot(index='y_coordinate', columns='x_coordinate', values=field)
    contour = ax.contourf(data.columns, data.index, data.values, cmap='viridis')
    ax.set_title(title, fontsize=14)
    ax.set_xlabel('$x$ (m)', fontsize=12)
    ax.set_ylabel('$y$ (m)', fontsize=12)
    fig.colorbar(contour, ax=ax, format='%.2e')
    return ax


def plot_ellipse(ax, a, b):
    ellipse = Ellipse((0, 0), 2 * a, 2 * b, edgecolor='black', facecolor='orange', alpha=0.3,
                      linewidth=2, linestyle='dashed')
    ax.add_patch(ellipse)
    ax.set_aspect('equal')
    return ax


if __name__ == "__main__":
    file_name = 'example_216_kleopatra.csv'
    df = read_and_process_data(file_name)

    # Ellipsoid parameters for Kleopatra
    a = 135000.0  # Semi-major axis
    b = 58000.0  # Semi-minor axis

    fig, axs = plt.subplots(1, 2, figsize=(15, 7))

    axs[0] = plot_field(df, 'potential', axs[0], 'Potential Field (XY Plane)')
    axs[0] = plot_ellipse(axs[0], a, b)

    axs[1] = plot_field(df, 'acceleration_magnitude', axs[1], 'Acceleration Field (XY Plane)')
    axs[1] = plot_ellipse(axs[1], a, b)

    plt.tight_layout()
    plt.show()
