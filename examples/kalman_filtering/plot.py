# Import libraries at the top
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import patches
import numpy as np

# Constants related to loading data
START_INDEX = 0
END_INDEX = -1
FILTER_NAN = True

# Constants related to opacities
INVERT_OPACITY = True
MIN_OPACITY = 0.02
MAX_OPACITY = 0.2

# Constants related to blow up detection
N_POINTS = 100  # Number of points to plot before blow up

# Load data function
def load_data(path, start, end, filter_nan=False):
    df = pd.read_csv(path).iloc[start:end]
    if filter_nan:
        df = df.dropna()
    return df


def plot_stations(ax, stations):
    for name, (x, y) in stations.items():
        ax.plot(x, y, 'ko')
        ax.text(x, y, name, fontsize=12, ha='right')


def plot_earth(ax):
    earth = patches.Circle((0, 0), 6371, fill=False, color='green', linestyle=':')
    ax.add_patch(earth)


# Function for calculating covariance ellipse parameters
def calculate_covariance_parameters(df):
    parameters = []

    for i, row in df.iterrows():
        cov = [[row['cov_xx'], row['cov_xy']], [row['cov_xy'], row['cov_yy']]]
        eigenvalues, eigenvectors = np.linalg.eig(cov)
        width, height = 2 * np.sqrt(eigenvalues)
        angle = np.degrees(np.arctan2(*eigenvectors[:, 0][::-1]))
        parameters.append((row['estimate_x'], row['estimate_y'], width, height, angle))

    return parameters


# Function for plotting the covariance ellipses
def plot_covariance_ellipses(
        ax, df, parameters, start, i_blow_up=None):
    color_map = plt.cm.get_cmap('viridis')

    for i, (estimate_x, estimate_y, width, height, angle) in enumerate(parameters):
        # Only draw the ellipse if it's within the last N_POINTS before the blowup
        if i_blow_up is None or (i <= i_blow_up - N_POINTS):
            color = color_map(i / len(parameters))
            ellipse = patches.Ellipse((estimate_x, estimate_y), width, height, angle=angle, fill=True, color=color, alpha=0.1)
            ax.add_patch(ellipse)



# Functions for plotting
# Functions for plotting
def plot_trajectory(df, ax, color='black', linestyle='--', label='Estimated', highlight_last_point=False):
    ax.plot(df['actual_x'], df['actual_y'], label='Actual', color='blue')
    ax.plot(df['estimate_x'], df['estimate_y'], label=label, color=color, linestyle=linestyle)

    # Highlight the last point if required
    if highlight_last_point:
        ax.plot(df.iloc[-1]['estimate_x'], df.iloc[-1]['estimate_y'], 'x', color='red', markersize=10)


def blow_up_detection(df, start, blow_up_threshold=0.03):
    prev_x, prev_y = df.iloc[start]['estimate_x'], df.iloc[start]['estimate_y']
    blow_up_detected = False
    i_blow_up = 0

    for i, row in df.iterrows():
        if i > start:
            current_x, current_y = row['estimate_x'], row['estimate_y']
            prev_norm = np.linalg.norm([prev_x, prev_y])
            current_norm = np.linalg.norm([current_x, current_y])

            if abs((current_norm - prev_norm) / prev_norm) > blow_up_threshold:
                print(f"System blew up at iteration: {i}")
                blow_up_detected = True
                i_blow_up = i
                break

        prev_x, prev_y = row['estimate_x'], row['estimate_y']

    return blow_up_detected, i_blow_up

WINDOW_START = -N_POINTS  # Change this to control the start of your window
WINDOW_END = -1  # Change this to control the end of your window, None for complete end
N_POINTS = 200  # Number of points before the blow-up to display
def main():
    filenames = [
        'cmake-build-release/trajectory_float_set1.csv',
        'cmake-build-release/trajectory_double_set1.csv',
        'cmake-build-release/trajectory_long_double_set1.csv'
    ]
    figure, axs = plt.subplots(1, len(filenames), figsize=(15, 5))  # adjust the figure size if needed

    for i, filename in enumerate(filenames):
        ax = axs[i]  # select the i-th subplot
        print(f"Loading data from {filename}...")
        df = load_data(filename, START_INDEX, END_INDEX, FILTER_NAN)

        # Blow-up detection
        blow_up_detected, i_blow_up = blow_up_detection(df, df.first_valid_index())

        if blow_up_detected:
            # Reduce the dataframe to include only the N_POINTS leading up to blow-up
            df = df[max(0, i_blow_up - N_POINTS):i_blow_up-2]

            print(f"System blew up, keeping only the {N_POINTS} iterations leading up to the blow-up.")
        else:
            print("System did not blow up. No plot generated.")

        # Adjust the data to the specified window
        df = df[WINDOW_START:WINDOW_END]

        # Calculate covariance ellipse parameters for the reduced dataframe
        parameters = calculate_covariance_parameters(df)

        # Plot covariance ellipses
        plot_covariance_ellipses(ax, df, parameters, df.first_valid_index(), i_blow_up=None)

        # Plot trajectory with highlighted last point
        plot_trajectory(df, ax, color='red', linestyle='--', label='Updated Estimate', highlight_last_point=True)

        estimate_x, estimate_y, width, height, angle = parameters[-1]
        ax.plot(df.tail(1)['estimate_x'], df.tail(1)['estimate_y'], 'ro', markersize=5)
        ellipse = patches.Ellipse((df.tail(1)['estimate_x'], df.tail(1)['estimate_y']), width, height, angle=angle,
                                  fill=False, color='orange')
        ax.add_patch(ellipse)

        stations = {'ts1': [-6371, 0], 'ts2': [6371, 0], 'ts3': [0, 6371], 'ts4': [0, -6371]}
        plot_stations(ax, stations)
        plot_earth(ax)

        ax.axis('equal')
        ax.set_title(filename)  # set the title of the plot to the filename
        ax.legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
