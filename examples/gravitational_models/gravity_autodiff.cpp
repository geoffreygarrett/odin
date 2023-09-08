#include <Eigen/Dense>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/real.hpp>
#include <cmath>
#include <iostream>
#include <odin/models/gravitational/spherical.hpp>

using namespace autodiff;

template<typename Function, typename Scalar, int N>
auto compute_finite_difference_gradient(
        Function func, const Eigen::Matrix<Scalar, N, 1>& point, double h = 1e-5) {
    Eigen::Matrix<Scalar, N, 1> gradient;
    using VectorN = Eigen::Matrix<Scalar, N, 1>;
    for (int i = 0; i < N; ++i) {
        VectorN point_plus  = point;
        VectorN point_minus = point;
        point_plus[i] += h;
        point_minus[i] -= h;

        auto f_plus  = func(point_plus);
        auto f_minus = func(point_minus);

        gradient[i] = (f_plus - f_minus) / (2.0 * h);
    }
    return gradient;
}


template<typename Function, typename Scalar, int N>
auto compute_finite_difference_jacobian(
        Function func, const Eigen::Matrix<Scalar, N, 1>& point, double h = 1e-5) {
    Eigen::Matrix<Scalar, N, N> jacobian;
    using VectorN = Eigen::Matrix<Scalar, N, 1>;
    for (int i = 0; i < N; ++i) {
        VectorN point_plus  = point;
        VectorN point_minus = point;
        point_plus[i] += h;
        point_minus[i] -= h;

        auto f_plus  = func(point_plus);
        auto f_minus = func(point_minus);

        // Compute the i-th column of the Jacobian
        VectorN column  = (f_plus - f_minus) / (2.0 * h);
        jacobian.col(i) = column;
    }
    return jacobian;
}

int main(int argc, char** argv) {
    using scalar_type  = autodiff::real;
    using vector3_type = Eigen::Vector<scalar_type, 3>;
    auto my_point_mass = odin::gravity::PointMass<scalar_type>(1.0, vector3_type(0.0, 0.0, 0.0));

    vector3_type test_point(1.0, 0.0, 0.0);

    // Test partial derivatives of Potential w.r.t position
    auto autodiff_partial_potential    = my_point_mass.partial_potential_wrt_position(test_point);
    auto finite_diff_partial_potential = compute_finite_difference_gradient(
            [&](auto&& point) { return my_point_mass.potential(point); }, test_point);

    // Test partial derivatives of Acceleration w.r.t position
    auto autodiff_partial_accel    = my_point_mass.partial_acceleration_wrt_position(test_point);
    auto finite_diff_partial_accel = compute_finite_difference_jacobian(
            [&](auto& point) { return my_point_mass.acceleration(point); }, test_point);

    // Compare results
    double tolerance = 1e-6;

    std::cout << "Autodiff partial potential: " << autodiff_partial_potential.transpose()
              << std::endl;
    std::cout << "Finite difference partial potential: "
              << finite_diff_partial_potential.transpose() << std::endl;
    std::cout << "Autodiff partial acceleration: " << autodiff_partial_accel.transpose()
              << std::endl;
    std::cout << "Finite difference partial acceleration: " << finite_diff_partial_accel.transpose()
              << std::endl;

    if ((autodiff_partial_potential - finite_diff_partial_potential).norm() < tolerance) {
        std::cout << "Partial derivatives of Potential Test passed! Autodiff and finite difference "
                     "results are close enough."
                  << std::endl;
    } else {
        std::cout << "Partial derivatives of Potential Test failed! Results are not close enough."
                  << std::endl;
    }
    //
    if ((autodiff_partial_accel - finite_diff_partial_accel).norm() < tolerance) {
        std::cout << "Partial derivatives of Acceleration Test passed! Autodiff and finite "
                     "difference results are close enough."
                  << std::endl;
    } else {
        std::cout
                << "Partial derivatives of Acceleration Test failed! Results are not close enough."
                << std::endl;
    }

    return 0;
}
