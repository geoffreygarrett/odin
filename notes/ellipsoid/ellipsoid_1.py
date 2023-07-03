import numpy as np
import matplotlib.pyplot as plt


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
        np.linspace(
            -2*(c**2), np.max([np.abs(x * a), np.abs(y * b), np.abs(z * c)]), 1000
        ),
        # np.linspace(-(c**2), -(b**2), 1000),
        # np.linspace(-(b**2), -(a**2), 1000),
        # np.linspace(-(a**2), -np.max([x * a, y * b, z * c]) ** 2.5, 1000),
    ]
    C_segment_values = []
    for segment in segments:
        C_values = C(x, y, z, segment, a, b, c)
        C_segment_values.append(C_values)
    return segments, C_segment_values


def plot_segments(kappa_segments, C_value_segments):
    """Function to plot the segments with different colors"""
    colors = ["blue", "green", "red"]
    labels = ["Ellipsoidal Confocal Family", "Hyperboloid 1", "Hyperboloid 2"]
    for kappa_values, C_values, color, label in zip(
        kappa_segments, C_value_segments, colors, labels
    ):
        plt.plot(kappa_values, C_values, color=color, label=label)

        # Calculate crossing points of the x-axis and plot the roots
        crossing_indices = np.where(np.diff(np.sign(C_values)))[0]
        roots = kappa_values[crossing_indices]
        print(roots)
        # for roots with two values, remove the first
        if len(roots) > 1:
            roots = np.delete(roots, np.where(roots == roots[0]))
        # for
        plt.plot(roots, np.zeros(len(roots)), "ro")
        # for root in roots:
        #     plt.annotate(
        #         f"{root:.2f}",
        #         (root, 0),
        #         textcoords="offset points",
        #         xytext=(0, 10),
        #         ha="center",
        #     )

    plt.title("Plot of C(kappa) vs kappa")
    plt.xlabel("kappa")
    plt.ylabel("C(kappa)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()


def main():
    # a: 300 b: 200 c: 100 x: -270 y: -10 z: 0 lower_bound: -10000 upper_bound: 6.561e+09
    # Set your ellipsoid constants
    a, b, c = 300, 200, 100

    # Set your point coordinates
    # x, y, z = -10, -10, -10
    x, y, z = -100, 0, 0
    # x, y, z = -300, 0, 0
    # x, y, z = 300,300,300

    kappa_segments, C_value_segments = calculate_segments(x, y, z, a, b, c)
    plt.figure(figsize=(10, 6))

    plot_segments(kappa_segments, C_value_segments)
    # verticle lines at each axis square
    plt.axvline(x=-(c**2), color="black", linestyle="--")
    # plt.axvline(x=-(b**2), color="black", linestyle="--")
    # plt.axvline(x=-(a**2), color="black", linestyle="--")
    plt.show()


if __name__ == "__main__":
    main()
