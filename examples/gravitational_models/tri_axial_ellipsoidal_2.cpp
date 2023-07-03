/**
* @file example_216_ellipsoid.cpp
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
#include <highfive/H5Easy.hpp>
#include <highfive/H5File.hpp>
#include <iomanip>// Needed for setprecision
#include <iostream>
#include <odin/logging.hpp>                       // Custom logging utilities
#include <odin/models/gravitational/spherical.hpp>// Custom class for modelling a sphere
//#include <odin/models/gravitational/spherical_harmonic.hpp>// Custom class for modelling a sphere
//#include <odin/models/gravitational/point_mass.hpp>            // Custom class for modelling a sphere
#include <odin/models/gravitational/ellipsoidal.hpp>// Custom class for modelling an ellipsoid
#include <tbb/tbb.h>
#include <vector>
// Define macro for changing tbb::concurrent_vector to std::vector
#define CONVERT_TO_STD_VECTOR(name) std::vector<decltype(name)::value_type>(name.begin(), name.end())
using Scalar         = double;
constexpr size_t dim = 3;
using Vector         = Eigen::Vector<Scalar, dim>;

struct GridPoint {
    Vector position;
    Scalar potential;
    Vector acceleration;
};

auto save_to_hd5 = [](
                           const tbb::concurrent_vector<GridPoint> &grid_points,
                           const std::string                       &filename,
                           double a, double b, double c,
                           double rho, double mu,
                           double step_size, double grid_range) {
    ODIN_LOG_INFO << "Writing to " << filename;

    // Log and write to HDF5
    HighFive::File file(filename, HighFive::File::Overwrite | HighFive::File::Create);

    std::vector<Vector> positions_std, accelerations_std;
    std::vector<Scalar> potentials_std;

    // Copy into std::vector
    for (auto &grid_point: grid_points) {
        positions_std.emplace_back(grid_point.position);
        potentials_std.emplace_back(grid_point.potential);
        accelerations_std.emplace_back(grid_point.acceleration);
    }

    // Dump all grid points data to the file
    H5Easy::dump(file, "/grid_points/positions", positions_std, H5Easy::DumpMode::Overwrite);
    H5Easy::dump(file, "/grid_points/potentials", potentials_std, H5Easy::DumpMode::Overwrite);
    H5Easy::dump(file, "/grid_points/accelerations", accelerations_std, H5Easy::DumpMode::Overwrite);

    // Assign attributes
    H5Easy::dumpAttribute(file, "/grid_points/positions", "a", a);
    H5Easy::dumpAttribute(file, "/grid_points/positions", "b", b);
    H5Easy::dumpAttribute(file, "/grid_points/positions", "c", c);
    H5Easy::dumpAttribute(file, "/grid_points/positions", "rho", rho);
    H5Easy::dumpAttribute(file, "/grid_points/positions", "mu", mu);
    H5Easy::dumpAttribute(file, "/grid_points/positions", "step_size", step_size);
    H5Easy::dumpAttribute(file, "/grid_points/positions", "grid_range", grid_range);

    ODIN_LOG_INFO << "Finished writing the potential field to " << filename;
};


// Jacbobi integral


int main() {
    INIT_ODIN_LOGGING("test_tri_axial_ellipsoid", "./log/test_tri_axial_ellipsoid.log");

    const double a                         = 300;
    const double b                         = 200;
    const double c                         = 100;
    const double gram_per_cubic_centimeter = 2.8;
    const double rho                       = gram_per_cubic_centimeter * 1000;  // Density in kg/m³
    const double G                         = 6.67408 * 1e-11;                   // Gravitational constant in m³/(kg·s²)
    const double mu                        = G * a * b * c * M_PI * 4 / 3 * rho;// Gravitational parameter in m³/s²

    const double step_size = 5.0, grid_range = 1200.0;
    const int    num_steps = static_cast<int>((2 * grid_range) / step_size) + 1;


    // Create the point mass model
    // Create a vector of grid points
    // translate a,b,c to the equivalent spherical radius
    double a_spherical = std::pow(a * b * c, 1.0 / 3.0);
    double b_spherical = a_spherical;
    double c_spherical = a_spherical;
    double r_spherical = a_spherical;

    tbb::concurrent_vector<GridPoint> grid_points_point_mass;
    grid_points_point_mass.reserve(num_steps * num_steps);
    HomogeneousSphere<Scalar> spherical(mu, r_spherical);
    ODIN_LOG_INFO << "Computed the potential field on a grid of " << grid_points_point_mass.size() << " points";
    ODIN_LOG_INFO << "Computing the potential field for a point mass at the origin";

#define OMIT_REGION 0.0
#define DEG_2_RAD (M_PI / 180.0)
#define ROTATION_RATE (0.0 * DEG_2_RAD)

    tbb::parallel_for(
            tbb::blocked_range2d<int>(0, num_steps, 0, num_steps),
            [&](const tbb::blocked_range2d<int> &r) {
                // Create a new Ellipsoid for each thread ( the solvers are not thread-safe )
                thread_local auto tl_spherical = spherical.thread_local_copy();

                for (int i = r.rows().begin(); i != r.rows().end(); ++i) {
                    for (int j = r.cols().begin(); j != r.cols().end(); ++j) {
                        double x = (i - num_steps / 2.0) * step_size;
                        double y = (j - num_steps / 2.0) * step_size;
                        double z = 0.0;
                        if ((x * x / (a_spherical * a_spherical) + y * y / (b_spherical * b_spherical) + z * z / (c_spherical * c_spherical) > OMIT_REGION)) {
                            Eigen::Vector<double, 3> position(x, y, z);
                            double                   potential = tl_spherical.potential(position);
                            potential += 0.5 * ROTATION_RATE * ROTATION_RATE * (x * x + y * y);
                            Eigen::Vector<double, 3> acceleration = tl_spherical.acceleration(position);
                            // Usually with vectors, emplace_back is NOT thread-safe, and we would need
                            // to use a mutex lock. However, tbb::concurrent_vector
                            // is thread-safe and we can use it instead without the added hassle.
                            grid_points_point_mass.emplace_back(GridPoint{position, potential, acceleration});
                        }
                    }
                }
            });

    // translate a,b,c to the average, so this is spherical
    save_to_hd5(grid_points_point_mass, "point_mass.h5", a_spherical, b_spherical, c_spherical, rho, mu, step_size, grid_range);

    // Create a vector of grid points
    // Create the ellipsoidal model
    tbb::concurrent_vector<GridPoint> grid_points_ellipsoid;
    grid_points_ellipsoid.reserve(num_steps * num_steps);
    TriAxialEllipsoid<Scalar> ellipsoid(a, b, c, mu);
    // TODO: Root finding fails for ellipsoid with 3 identical axes. I need to change some assumptions in the routine.

    ODIN_LOG_INFO << "Computing the potential field on a grid of " << num_steps << " x " << num_steps << " points";
    tbb::parallel_for(
            tbb::blocked_range2d<int>(0, num_steps, 0, num_steps),
            [&](const tbb::blocked_range2d<int> &r) {
                // Create a new Ellipsoid for each thread ( the solvers are not thread-safe )
                thread_local auto tl_ellipsoid = ellipsoid.thread_local_copy();

                for (int i = r.rows().begin(); i != r.rows().end(); ++i) {
                    for (int j = r.cols().begin(); j != r.cols().end(); ++j) {
                        double x = (i - num_steps / 2.0) * step_size;
                        double y = (j - num_steps / 2.0) * step_size;
                        double z = 0.0;
                        if ((x * x / (a * a) + y * y / (b * b) + z * z / (c * c) > OMIT_REGION)) {
                            Eigen::Vector<double, 3> position(x, y, z);
                            double                   potential    = tl_ellipsoid.potential(position);
                            potential += 0.5 * ROTATION_RATE * ROTATION_RATE * (x * x + y * y);
                            Eigen::Vector<double, 3> acceleration = tl_ellipsoid.acceleration(position);
                            // Usually vectors' emplace_back is NOT thread-safe and we would need
                            // to use a mutex to protect the vector. However, tbb::concurrent_vector
                            // is thread-safe and we can use it instead without the added hassle.
                            grid_points_ellipsoid.emplace_back(GridPoint{position, potential, acceleration});
                        }
                    }
                }
            });

    save_to_hd5(grid_points_ellipsoid, "ellipsoid.h5", a, b, c, rho, mu, step_size, grid_range);


    return 0;
}
