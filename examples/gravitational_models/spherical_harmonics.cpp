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
    double potential;
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

//
//template<typename Scalar, int Degree = -1, int Order = -1>
//struct GeopotentialCoefficientTransformations {
//    using Clm = std::vector<std::vector<Scalar>>;
//    using Slm = std::vector<std::vector<Scalar>>;
//
//    static void add_dimensions_in_place(Scalar mu, Scalar R, Clm &C_lm, Slm &S_lm) {
//        for (int l = 0; l < C_lm.size(); l++) {
//            Scalar scale_factor = -mu * std::pow(R, l);
//            for (int m = 0; m < C_lm[l].size(); m++) {
//                C_lm[l][m] *= scale_factor;
//                S_lm[l][m] *= scale_factor;
//            }
//        }
//    }

template<typename Scalar, int Degree = 8, int Order = 8>
SphericalHarmonics<Scalar, Degree, Order> create_spherical_harmonics_model(Scalar mu, Scalar radius) {
    ODIN_LOG_INFO << "Creating the spherical harmonics model";

    // Define and Initialize your coefficients
    typename SphericalHarmonics<Scalar, Degree, Order>::Clm Cnm_values{};
    typename SphericalHarmonics<Scalar, Degree, Order>::Slm Snm_values{};

    //    const Scalar scale = - mu * pow(radius, 2) / 3.0;

    // Zonal coefficients
    Cnm_values[0][0] = -1.0;
    Snm_values[0][0] = 0.0;
    Cnm_values[2][0] = -0.1082635854e-2;
    Cnm_values[3][0] = 0.2532435346e-5;
    Cnm_values[4][0] = 0.1619331205e-5;
    Cnm_values[5][0] = 0.2277161016e-6;
    Cnm_values[6][0] = -0.5396484906e-6;
    Cnm_values[7][0] = 0.3513684422e-6;
    Cnm_values[8][0] = 0.2025187152e-6;
    //
    //    // Tesseral coefficients
    Cnm_values[2][1] = -0.3504890360e-9;
    Cnm_values[2][2] = 0.1574536043e-5;
    Cnm_values[3][1] = 0.2192798802e-5;
    Cnm_values[3][2] = 0.3090160446e-6;
    Cnm_values[3][3] = 0.1005588574e-6;
    Cnm_values[4][1] = -0.5087253036e-6;
    Cnm_values[4][2] = 0.7841223074e-7;
    Cnm_values[4][3] = 0.5921574319e-7;
    Cnm_values[4][4] = -0.3982395740e-8;
    //
    Snm_values[2][1] = 0.1635406077e-8;
    Snm_values[2][2] = -0.9038680729e-6;
    Snm_values[3][1] = 0.2680118938e-6;
    Snm_values[3][2] = -0.2114023978e-6;
    Snm_values[3][3] = 0.1972013239e-6;
    Snm_values[4][1] = -0.4494599352e-6;
    Snm_values[4][2] = 0.1481554569e-6;
    Snm_values[4][3] = -0.1201129183e-7;
    Snm_values[4][4] = 0.6525605810e-8;

    GeopotentialCoefficientTransformations<Scalar, Degree, Order>::add_dimensions_in_place(
            //            mu,
            398600.4415e9,
            6378.1363e3,
            Cnm_values, Snm_values);

    // Create an instance of SphericalHarmonics and assign the coefficients
    SphericalHarmonics<Scalar, Degree, Order> model(mu, radius, Cnm_values, Snm_values);
    //    SphericalHarmonics<Scalar, 8, 8> point_mass(mu, 8, a_spherical, spherical_harmonic);


    ODIN_LOG_INFO << "Finished creating the spherical harmonics model";

    return model;
};


// Standard factorial function
unsigned long long factorial(int n) {
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}


template<typename Scalar, int Degree, int Order = Degree>
SphericalHarmonics<Scalar, Degree, Order> create_spherical_harmonics_model(Scalar mu, Scalar R, Scalar a, Scalar b, Scalar c) {
    LOG_INFO << "Creating the spherical harmonics model";

    // Define and Initialize your coefficients
    typename SphericalHarmonics<Scalar, Degree, Order>::Clm Cnm_values{};
    typename SphericalHarmonics<Scalar, Degree, Order>::Slm Snm_values{};

    // calculate coefficients for C_{2l,2m}
    for (int l = 0; l <= Degree / 2; ++l) {
        for (int m = 0; m <= Order / 2; ++m) {
            Scalar sum_k = 0;
            for (int k = 0; k <= l - m; ++k) {
                Scalar sum_s = 0;
                for (int s = 0; s <= m; ++s) {
                    Scalar sum_p = 0;
                    for (int p = 0; p <= k; ++p) {
                        Scalar sum_q = 0;
                        for (int q = 0; q <= p; ++q) {
                            //                            sum_q += factorial(2 * s + 2 * p - 2 * q) * factorial(2 * m - 2 * s + 2 * q) / (factorial(q) * factorial(p - q) * factorial(s + p - q) * factorial(m - s + q)) * pow(a, 2 * (m - s + q)) * pow(b, 2 * (s + p - q)) * pow(c, 2 * (l - m - p));
                            sum_q += factorial(2 * s + 2 * p - 2 * q) * factorial(2 * m - 2 * s + 2 * q) / (factorial(q) * factorial(p - 1) * factorial(s + p - q) * factorial(m - s + q)) * pow(a, 2 * (m - s + q)) * pow(b, 2 * (s + p - q)) * pow(c, 2 * (l - m - p));
                        }
                        sum_p += factorial(2 * l - 2 * m - 2 * p) / (factorial(l - m - p) * factorial(k - p)) * sum_q;
                    }
                    sum_s += pow(-1, s) / (factorial(2 * s) * factorial(2 * m - 2 * s)) * sum_p;
                }
                sum_k += pow(-1, k) * factorial(4 * l - 2 * k) / (factorial(2 * l - k) * factorial(2 * l - 2 * m - 2 * k)) * sum_s;
            }
            Cnm_values[2 * l][2 * m] = pow(3.0 / (R * R * R * R), l) * factorial(l) * factorial(2 * m) / (pow(2, 2 * l) * (2 * l + 3) * factorial(2 * l + 1)) * (2 - (m == 0 ? 1 : 0)) * factorial(2 * l - 2 * m) / factorial(2 * l + 2 * m) * sum_k;
        }
    }

    SphericalHarmonics<Scalar, Degree, Order> model(mu, R, Cnm_values, Snm_values);
    //    model.set_coefficients(Cnm_values, Snm_values);

    LOG_INFO << "Finished creating the spherical harmonics model";

    return model;
};


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
    // translate a,b,c to the average, so this is spherical
    double a_spherical = pow(a * b * c, 1.0 / 3.0);
    double b_spherical = a_spherical;
    double c_spherical = a_spherical;

    //    auto spherical_harmonic = create_spherical_harmonics_model<Scalar, 10, 10>(mu, a_spherical, a, b, c);
    auto spherical_harmonic = create_spherical_harmonics_model(mu, a_spherical);

    tbb::concurrent_vector<GridPoint> grid_points_point_mass;
    grid_points_point_mass.reserve(num_steps * num_steps);

    //    SphericalHarmonics<Scalar, 8, 8> point_mass(mu, 8, a_spherical, spherical_harmonic);

    //    PointMass<Scalar> point_mass(mu);//Vector::Zero(),
    //    HollowSphere<Scalar> point_mass(mu, a_spherical);
    //    HomogeneousSphere<Scalar> point_mass(mu, a_spherical);
    //    IsochroneSphere<Scalar> point_mass(mu, a);

    ODIN_LOG_INFO << "Computed the potential field on a grid of " << grid_points_point_mass.size() << " points";
    ODIN_LOG_INFO << "Computing the potential field for a point mass at the origin";


#define OMIT_REGION 0.0
    tbb::parallel_for(
            tbb::blocked_range2d<int>(0, num_steps, 0, num_steps),
            [&](const tbb::blocked_range2d<int> &r) {
                // Create a new Ellipsoid for each thread ( the solvers are not thread-safe )
                thread_local auto tl_point_mass = spherical_harmonic.thread_local_copy();

                for (int i = r.rows().begin(); i != r.rows().end(); ++i) {
                    for (int j = r.cols().begin(); j != r.cols().end(); ++j) {
                        double x = (i - num_steps / 2.0) * step_size;
                        double y = (j - num_steps / 2.0) * step_size;
                        double z = 0.0;
                        if ((x * x / (a_spherical * a_spherical) + y * y / (b_spherical * b_spherical) + z * z / (c_spherical * c_spherical) > OMIT_REGION)) {
                            Eigen::Vector<double, 3> position(x, y, z);
                            double                   potential    = spherical_harmonic.potential(position);
                            Eigen::Vector<double, 3> acceleration = spherical_harmonic.acceleration(position);
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
