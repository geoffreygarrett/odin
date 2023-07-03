
#include <Eigen/Core>
#include <cmath>
#include <gtest/gtest.h>
#include <odin/core/numerics/integrators/runge_kutta.hpp>
#include <tuple>
// Function y′ = 1 + y^2
//auto func_y_prime = [](double y, double) {
//    return y * y + 1;
//};

auto func_y_prime = [](const Eigen::Vector<double, 1> &y, double) {
    Eigen::Vector<double, 1> result;
    result << y(0) * y(0) + 1;
    return result;
};

//std::function<Eigen::Vector<double, 1>(Eigen::Vector<double, 1> &, double)> func_y_prime_eigen = func_y_prime;

//auto func_y_prime = [](const Eigen::Vector<double, 1> &y, double) {
//    return y * y + 1;
//};

//std::function<Eigen::Vector<double, 1>(const Eigen::Vector<double, 1> &, double)> func_y_prime_eigen = func_y_prime;


/**
 * @brief Test data for the RK45 method.
 *
 * The test data comes from the following source:
 * https://maths.cnam.fr/IMG/pdf/RungeKuttaFehlbergProof.pdf (page 499, table 9.11)
 * The table is titled "RKF45 Solution to y′ = 1 + y^2 , y(0) = 0"
 *
 * The test data consists of tuples:
 *   k: step number
 *   t_k: time at step k
 *   y_k: RK45 approximation
 *   y_true: True solution, y(t k ) = tan(t k )
 *   error: y(t k ) − yk
 */
std::vector<std::tuple<int, double, double, double, double>> test_cases = {
  //k   t_k     y_k        y_true     error
        { 0,  0.0, 0.0000000, 0.0000000,  0.0000000},
        { 1,  0.2, 0.2027100, 0.2027100,  0.0000000},
        { 2,  0.4, 0.4227933, 0.4227931, -0.0000002},
        { 3,  0.6, 0.6841376, 0.6841368, -0.0000008},
        { 4,  0.8, 1.0296434, 1.0296386, -0.0000048},
        { 5,  1.0, 1.5574398, 1.5774077, -0.0000321},
        { 6,  1.1, 1.9648085, 1.9647597, -0.0000488},
        { 7,  1.2, 2.5722408, 2.5721516, -0.0000892},
        { 8,  1.3, 3.6023295, 3.6021024, -0.0002271},
        { 9, 1.35, 4.4555714, 4.4552218, -0.0003496},
        {10,  1.4, 5.7985045, 5.7978837, -0.0006208}
};

TEST(RK45Test, FunctionYPrimeTest) {
    // Initialize RK45 with the function
    RungeKuttaIntegrator<
            Eigen::Vector<double, 1>,
            decltype(func_y_prime),
            butcher_tableau::Fehlberg<double>>
            integrator(func_y_prime);

    Eigen::Vector<double, 1> y{0.0};    // initial condition
    const double             dt  = 0.1; // priori step size
    const double             tol = 2e-5;// tolerance
    //    y = 0;// y(0) = 0

    // Iterate over test data
    for (const auto &test_case: test_cases) {
        auto [k, tk, yk, y_true, error] = test_case;

        double current_time = 0.0;

        // Perform integration steps until reaching time tk
        while (current_time < tk) {
            std::pair(y, dt) = integrator.step_with_dt(y, 0.001);// Note: step size may need adjusting according to your implementation

            current_time += 0.001;
        }

        ASSERT_NEAR(y[0], yk, 1e-6);// Check if the RK45 approximation is correct within tolerance
    }
}