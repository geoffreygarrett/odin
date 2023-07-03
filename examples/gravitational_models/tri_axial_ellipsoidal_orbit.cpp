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

// Update your SpacecraftTrajectory struct
template<typename Scalar, size_t Dim>
struct SpacecraftTrajectory {
    std::vector<Eigen::Vector<Scalar, Dim>> positions;
    std::vector<Eigen::Vector<Scalar, Dim>> velocities;
    std::vector<Eigen::Vector<Scalar, Dim>> accelerations;
    std::vector<Scalar>                     times;
};

struct GridPoint {
    Vector position;
    Scalar potential;
    Vector acceleration;
};

auto save_gravitation_to_hd5 = [](
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
//HIGHFIVE_REGISTER_TYPE(SpacecraftTrajectory, SpacecraftTrajectory)

auto save_to_hd5 = [](
                           const std::vector<SpacecraftTrajectory<Scalar, 3>> &trajectories,
                           const std::string                                  &filename,
                           const double                                        a,
                           const double                                        b,
                           const double                                        c,
                           const double                                        rho,
                           const double                                        mu,
                           const double                                        step_size,
                           const double                                        sim_time) {
    ODIN_LOG_INFO << "Writing to " << filename;

    // Create an HDF5 file
    HighFive::File file(filename, HighFive::File::Overwrite | HighFive::File::Create);

    // Create a new group in the file to contain the trajectories
    HighFive::Group group = file.createGroup("/trajectories");

    // Iterate over all the trajectories and dump each to the file
    for (size_t i = 0; i < trajectories.size(); i++) {
        // Convert i to string for creating distinct dataset names
        std::string dataset_name = "trajectory_" + std::to_string(i);

        // Create a group for each trajectory
        HighFive::Group traj_group = group.createGroup(dataset_name);

        // Save each part of the trajectory as a separate dataset in the group
        traj_group.createDataSet("positions", trajectories[i].positions);
        traj_group.createDataSet("velocities", trajectories[i].velocities);
        traj_group.createDataSet("times", trajectories[i].times);
    }

    // Assign attributes to the group
    group.createAttribute<Scalar>("a", a);
    group.createAttribute<Scalar>("b", b);
    group.createAttribute<Scalar>("c", c);
    group.createAttribute<Scalar>("rho", rho);
    group.createAttribute<Scalar>("mu", mu);
    group.createAttribute<Scalar>("step_size", step_size);
    group.createAttribute<Scalar>("sim_time", sim_time);

    //    H5Easy::dumpAttribute(file, "/trajectories", "a", a);
    //    H5Easy::dumpAttribute(file, "/trajectories", "b", b);
    //    H5Easy::dumpAttribute(file, "/trajectories", "c", c);
    //    H5Easy::dumpAttribute(file, "/trajectories", "rho", rho);
    //    H5Easy::dumpAttribute(file, "/trajectories", "mu", mu);
    //    H5Easy::dumpAttribute(file, "/trajectories", "step_size", step_size);
    //    H5Easy::dumpAttribute(file, "/trajectories", "sim_time", sim_time);

    ODIN_LOG_INFO << "Finished writing the potential field to " << filename;
};
#include "odin/domain/astrodynamics.hpp"

#define DEG_TO_RAD (M_PI / 180.0)
#define RAD_TO_DEG (180.0 / M_PI)


// Compute acceleration in inertial frame given time, position and rotational speed
Vector compute_inertial_acceleration(Scalar time, Vector position, TriAxialEllipsoid<Scalar> &ellipsoid, Scalar rate_of_rotation) {
    // Compute the model-fixed acceleration
    const Vector model_center_position    = Eigen::Vector3d::Zero();
    const auto   model_rotation_matrix    = Eigen::AngleAxis<Scalar>(rate_of_rotation * time, Vector::UnitZ()).toRotationMatrix();
    const Vector model_fixed_position     = model_rotation_matrix.transpose() * (position - model_center_position);
    const Vector model_fixed_acceleration = ellipsoid.acceleration(model_fixed_position);

    // Transform the acceleration to the inertial frame
    const Vector inertial_acceleration = model_rotation_matrix * model_fixed_acceleration;

    return inertial_acceleration;
}

int main() {
    INIT_ODIN_LOGGING("test_tri_axial_ellipsoid", "./log/test_tri_axial_ellipsoid.log");

    const double a                         = 300;
    const double b                         = 200;
    const double c                         = 100;
    const double gram_per_cubic_centimeter = 2.8;
    const double rho                       = gram_per_cubic_centimeter * 1000;  // Density in kg/m³
    const double G                         = 6.67408 * 1e-11;                   // Gravitational constant in m³/(kg·s²)
    const double mu                        = G * a * b * c * M_PI * 4 / 3 * rho;// Gravitational parameter in m³/s²
    // Create the point mass model
    // Create a vector of grid points
    // translate a,b,c to the average, so this is spherical
    //    double                    a_spherical = (a + b + c) / 3.0;
    //    double                    b_spherical = (a + b + c) / 3.0;
    //    double                    c_spherical = (a + b + c) / 3.0;
    //    HomogeneousSphere<Scalar> point_mass(mu, a_spherical);
    TriAxialEllipsoid<Scalar> ellipsoid(a, b, c, mu);

    // Create initial conditions for the sets of trajectories to be simulated. We will define it in Keplerian elements and convert it to Cartesian coordinates
    std::vector<SpacecraftTrajectory<Scalar, 3>> trajectories;
    const Scalar                                 t_start = 0.0;
    const Scalar                                 t_end   = 2 * M_PI * sqrt(600.0 * 600.0 * 600.0 / mu);
    const Scalar                                 dt      = 5.0;
    const int                                    steps   = static_cast<int>(t_end / dt) + 1;
    trajectories.reserve(steps);
    for (int i = 0; i < steps; ++i) {
        SpacecraftTrajectory<Scalar, 3> trajectory;
        trajectory.positions.reserve(steps);
        trajectory.velocities.reserve(steps);
        trajectory.accelerations.reserve(steps);
        trajectory.times.reserve(steps);

        // Define the initial conditions in Keplerian elements
        const double sma   = 600.0 + 10.0 * i / 10.0;
        const double e     = 0.0;
        const double inc   = 0.0;
        const double Omega = 0.0;
        const double omega = 0.0;
        const double nu    = 0.0;
        const double p     = sma * (1.0 - e * e);

        // Convert the initial conditions to Cartesian coordinates
        auto [position, velocity] = coe2rv(mu, p, e, inc * DEG_TO_RAD, Omega * DEG_TO_RAD, omega * DEG_TO_RAD, nu * DEG_TO_RAD);
        const double time         = t_start;


        // Add the initial conditions to the trajectory
        trajectory.positions.push_back(position);
        trajectory.velocities.push_back(velocity);
        trajectory.times.push_back(time);

        // Add the trajectory to the vector of trajectories
        trajectories.push_back(trajectory);
    }

    // rate of rotation is 1 degree per second
    const double rate_of_rotation = 0.0000 * DEG_TO_RAD;

    // simulate a single spacecraft trajectory with rk4 for verification
    for (int i = 1; i < steps; ++i) {
        // Get the current position and velocity
        const Vector &position = trajectories[0].positions.back();
        const Vector &velocity = trajectories[0].velocities.back();
        const Scalar  time     = trajectories[0].times.back();

        // Compute the acceleration
        //        const Vector model_center_position    = Eigen::Vector3d::Zero();
        //        const auto   model_rotation_matrix    = Eigen::AngleAxis<Scalar>(rate_of_rotation * time, Vector::UnitZ()).toRotationMatrix();
        //        const Vector model_fixed_position     = model_rotation_matrix.transpose() * (position - model_center_position);
        //        const Vector model_fixed_acceleration = ellipsoid.acceleration(model_fixed_position);
        //        const Vector inertial_acceleration    = model_rotation_matrix * model_fixed_acceleration;
        //
        //        std::cout << "t: " << time << " rot: " << rate_of_rotation * time * RAD_TO_DEG << std::endl;


        // RK4 integration step
        // Compute k1
        //        Vector compute_inertial_acceleration(Scalar time, Vector position, TriAxialEllipsoid<Scalar> &ellipsoid, Scalar rate_of_rotation) {
        Vector k1_v = compute_inertial_acceleration(time, position, ellipsoid, rate_of_rotation);
        Vector k1_r = velocity;

        // Compute k2
        Vector r_temp = position + 0.5 * dt * k1_r;
        Vector v_temp = velocity + 0.5 * dt * k1_v;
        Vector k2_v   = compute_inertial_acceleration(time + 0.5 * dt, r_temp, ellipsoid, rate_of_rotation);
        //        Vector k2_v   = ellipsoid.acceleration(r_temp);
        Vector k2_r = v_temp;

        // Compute k3
        r_temp = position + 0.5 * dt * k2_r;
        v_temp = velocity + 0.5 * dt * k2_v;
        //        Vector k3_v = ellipsoid.acceleration(r_temp);
        Vector k3_v = compute_inertial_acceleration(time + 0.5 * dt, r_temp, ellipsoid, rate_of_rotation);
        Vector k3_r = v_temp;

        // Compute k4
        r_temp = position + dt * k3_r;
        v_temp = velocity + dt * k3_v;
        //        Vector k4_v = ellipsoid.acceleration(r_temp);
        Vector k4_v = compute_inertial_acceleration(time + dt, r_temp, ellipsoid, rate_of_rotation);
        Vector k4_r = v_temp;

        // Compute the new position and velocity
        Vector new_position = position + (dt / 6.0) * (k1_r + 2.0 * k2_r + 2.0 * k3_r + k4_r);
        Vector new_velocity = velocity + (dt / 6.0) * (k1_v + 2.0 * k2_v + 2.0 * k3_v + k4_v);

        // Add the new position and velocity to the trajectory
        trajectories[0].positions.push_back(new_position);
        trajectories[0].velocities.push_back(new_velocity);
        trajectories[0].times.push_back(i * dt);
    }

    //    auto save_to_hd5 = [](
    //                           const std::vector<SpacecraftTrajectory<Scalar, 3>> &trajectories,
    //                           const std::string                                  &filename,
    //                           const double                                        a,
    //                           const double                                        b,
    //                           const double                                        c,
    //                           const double                                        rho,
    //                           const double                                        mu,
    //                           const double                                        dt,
    //                           const double                                        sim_time) {
    save_to_hd5(trajectories, "test.h5", a, b, c, rho, mu, dt, steps);


    // The 'blocked_range' object partitions the problem size, here trajectories[0].positions.rows() - 1
    //    tbb::parallel_for(
    //            tbb::blocked_range<int>(1, 1000),
    //            [&](const tbb::blocked_range<int> &range) {
    //                for (int i = range.begin(); i != range.end(); ++i) {
    //                    // Get the current position and velocity
    //                    const Vector &position = trajectories[0].positions.back();
    //                    const Vector &velocity = trajectories[0].velocities.back();
    //
    //                    // Compute the acceleration
    //                    const Vector acceleration = ellipsoid.acceleration(position);
    //
    //                    // RK4 integration step
    //                    // Compute k1
    //                    Vector k1_v = acceleration;
    //                    Vector k1_r = velocity;
    //
    //                    // Compute k2
    //                    Vector r_temp = position + 0.5 * dt * k1_r;
    //                    Vector v_temp = velocity + 0.5 * dt * k1_v;
    //                    Vector k2_v   = ellipsoid.acceleration(r_temp);
    //                    Vector k2_r   = v_temp;
    //
    //                    // Compute k3
    //                    r_temp      = position + 0.5 * dt * k2_r;
    //                    v_temp      = velocity + 0.5 * dt * k2_v;
    //                    Vector k3_v = ellipsoid.acceleration(r_temp);
    //                    Vector k3_r = v_temp;
    //
    //                    // Compute k4
    //                    r_temp      = position + dt * k3_r;
    //                    v_temp      = velocity + dt * k3_v;
    //                    Vector k4_v = ellipsoid.acceleration(r_temp);
    //                    Vector k4_r = v_temp;
    //
    //                    // Compute the new position and velocity
    //                    Vector new_position = position + (dt / 6.0) * (k1_r + 2.0 * k2_r + 2.0 * k3_r + k4_r);
    //                    Vector new_velocity = velocity + (dt / 6.0) * (k1_v + 2.0 * k2_v + 2.0 * k3_v + k4_v);
    //
    //                    // Assign the new position and velocity to the trajectory
    //                    trajectories[0].positions.row(i)  = new_position;
    //                    trajectories[0].velocities.row(i) = new_velocity;
    //                    trajectories[0].times(i)          = i * dt;
    //                }
    //            });

    //        const double a                         = 300;
    //        const double b                         = 200;
    //        const double c                         = 100;
    //        const double gram_per_cubic_centimeter = 2.8;
    //        const double rho                       = gram_per_cubic_centimeter * 1000;  // Density in kg/m³
    //        const double G                         = 6.67408 * 1e-11;                   // Gravitational constant in m³/(kg·s²)
    //        const double mu                        = G * a * b * c * M_PI * 4 / 3 * rho;// Gravitational parameter in m³/s²

    const double step_size = 1.0, grid_range = 1000.0;
    const int    num_steps = static_cast<int>((2 * grid_range) / step_size) + 1;


    // Create a vector of grid points
    tbb::concurrent_vector<GridPoint> grid_points;
    grid_points.reserve(num_steps * num_steps);

#define OMIT_REGION 0.0

    // Create the ellipsoidal model
    tbb::concurrent_vector<GridPoint> grid_points_ellipsoid;
    grid_points_ellipsoid.reserve(num_steps * num_steps);
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
                            Eigen::Vector<double, 3> acceleration = tl_ellipsoid.acceleration(position);
                            // Usually vectors' emplace_back is NOT thread-safe and we would need
                            // to use a mutex to protect the vector. However, tbb::concurrent_vector
                            // is thread-safe and we can use it instead without the added hassle.
                            grid_points_ellipsoid.emplace_back(GridPoint{position, potential, acceleration});
                        }
                    }
                }
            });

    save_gravitation_to_hd5(grid_points_ellipsoid, "ellipsoid.h5", a, b, c, rho, mu, step_size, grid_range);


    return 0;
}
