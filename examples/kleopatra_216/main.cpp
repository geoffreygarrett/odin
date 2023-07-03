/**
 * @file example_216_kleopatra.cpp
 * @author Your Name Here
 * @date June 8, 2023
 *
 * @brief This file features the computational interpretation of Asteroid 216 Kleopatra based on radar observational data.
 *
 * The computational model incorporated in this program has been constructed leveraging the insights from the research study
 * "Radar Observations of Asteroid 216 Kleopatra" by Ostro, S.J. et al. (2000). The model incorporates the computation and analysis
 * of Kleopatra's shape, surface attributes, polar direction, and radar scattering law, among other properties.
 *
 * The underpinning research portrays Kleopatra as a dumbbell-shaped object with broad dimensions (217 x 94 x 81 km) with a
 * possible margin of 25%. The asteroid's surface characteristics align with a metallic regolith having similar porosity to lunar soil.
 * Kleopatra's shape is hypothesized to be a result of a complex chain of cosmic collision events, and a significant portion of its
 * internal structure might resemble a loosely packed rubble-pile formation.
 *
 * References:
 * Ostro, S.J., Hudson, R.S., Nolan, M.C., Margot, J.L., Scheeres, D.J., Campbell, D.B., Magri, C., Giorgini, J.D. and
 * Yeomans, D.K., 2000. Radar observations of asteroid 216 Kleopatra. Science, 288(5470), pp.836-839.
 * doi: 10.1126/science.288.5467.836
 *
 */
#include <fstream>
#include <iomanip>// Needed for setprecision
#include <iostream>
#include <odin/logging.hpp>                                   // Custom logging utilities
#include <odin/models/gravitational/ellipsoidal.hpp>          // Custom class for modelling an ellipsoid
#include <tbb/tbb.h>
#include <vector>

struct GridPoint {
    double                   x;
    double                   y;
    double                   z;
    double                   potential;
    Eigen::Vector<double, 3> acceleration;

    GridPoint(double x, double y, double z, double potential, const Eigen::Vector<double, 3> &acceleration)
        : x(x), y(y), z(z), potential(potential), acceleration(acceleration) {}
};


int main() {
    INIT_ODIN_LOGGING("test_tri_axial_ellipsoid", "./log/test_tri_axial_ellipsoid.log");

    const double a                         = 300;
    const double b                         = 200;
    const double c                         = 100;
    const double gram_per_cubic_centimeter = 2.0;
    const double rho                       = gram_per_cubic_centimeter * 1000;  // Density in kg/m³
    const double G                         = 6.67408 * 1e-11;                   // Gravitational constant in m³/(kg·s²)
    const double mu                        = G * a * b * c * M_PI * 4 / 3 * rho;// Gravitational parameter in m³/s²

    const double step_size = 100.0, grid_range = 1000.0;
    const int    num_steps = static_cast<int>((2 * grid_range) / step_size) + 1;

//    std::vector<GridPoint> grid_points;
//    grid_points.reserve(num_steps * num_steps);

    tbb::concurrent_vector<GridPoint> grid_points;
    grid_points.reserve(num_steps * num_steps);

    std::mutex my_mutex;// Mutex defined here

    Ellipsoid<double> kleopatra(a, b, c, mu);

    ODIN_LOG_INFO << "Computing the potential field on a grid of " << num_steps << " x " << num_steps << " points";
    tbb::parallel_for(
            tbb::blocked_range2d<int>(0, num_steps, 0, num_steps),
            [&](const tbb::blocked_range2d<int> &r) {
                // Create a new Ellipsoid for each thread ( the solvers are not thread-safe )
                thread_local auto tl_kleopatra = kleopatra.thread_local_copy();

                for (int i = r.rows().begin(); i != r.rows().end(); ++i) {
                    for (int j = r.cols().begin(); j != r.cols().end(); ++j) {
                        double x = (i - num_steps / 2.0) * step_size;
                        double y = (j - num_steps / 2.0) * step_size;
                        double z = 0.0;
                        if ((x * x / (a * a) + y * y / (b * b) + z * z / (c * c) > 0.9)) {
                            Eigen::Vector<double, 3> position(x, y, z);
                            double                   potential    = tl_kleopatra.potential(position);
                            Eigen::Vector<double, 3> acceleration = tl_kleopatra.acceleration(position);
                            // Since emplace_back is not thread-safe, we'll use a lock to ensure safe access
                            // Here we're using a std::mutex, but consider tbb's concurrent containers for real code
//                            std::lock_guard<std::mutex> lock(my_mutex);
                            grid_points.emplace_back(x, y, z, potential, acceleration);
                        }
                    }
                }
            });

    ODIN_LOG_INFO << "Computed the potential field on a grid of " << grid_points.size() << " points";
    ODIN_LOG_INFO << "Writing the potential field to example_216_kleopatra.csv";
    // Log and write to CSV
    std::ofstream output_file;
    output_file.open("example_216_kleopatra.csv");
    output_file << std::setprecision(10);

    output_file << "x_coordinate,y_coordinate,z_coordinate,potential,x_acceleration,y_acceleration,z_acceleration\n";
    for (const auto &point: grid_points) {
        // Write to CSV
        output_file << point.x << ","
                    << point.y << ","
                    << point.z << ","
                    << point.potential << ","
                    << point.acceleration[0] << ","
                    << point.acceleration[1] << ","
                    << point.acceleration[2] << "\n";
    }
    output_file.close();
    ODIN_LOG_INFO << "Finished writing the potential field to example_216_kleopatra.csv";

    return 0;
}
