#include <Eigen/Core>
#include <cmath>
#include <gtest/gtest.h>
#include <odin/core/numerics/integrators/runge_kutta.hpp>
#include <tuple>

using namespace butcher_tableau;

constexpr double t_final = 20.0;

template<typename Function, typename... Args>
using ReturnType = decltype(std::declval<Function>()(std::declval<Args>()...));


// Define a test fixture class.
template<class T>
class RungeKuttaIntegratorTest : public ::testing::Test {
public:
    Eigen::VectorXd initial;

    // Define some common setup, variables here.
    double dt        = 0.001;
    double t_final   = 20.0;
    double t_current = 0.0;

    template<typename Func>
    void TestIntegrator(Func func, double expected, int dim) {
        using IntegratorType = RungeKuttaIntegrator<Eigen::VectorXd, Func, T>;
        initial              = Eigen::VectorXd(dim);
        initial.setOnes();
        IntegratorType integrator(func);
        while (t_current < t_final) {
            auto x  = integrator.step(initial, dt);
            initial = x;
            t_current += dt;
            //            std::cout << t_current << ": " << x << std::endl;
            //            LOG(INFO) << t_current << ": " << x.transpose();
        }
        auto final = integrator.integrate(initial, dt, t_final);

        //        LOG(INFO) << "Final: " << final.transpose();
        EXPECT_NEAR(final[0], expected, 1e-6);
    }
};


// The types are now the fully defined integrator types
using IntegratorTypes = ::testing::Types<
        //        Euler<double>,   // TODO: THIS IS BROKEN xD, beyond its poor performance
        Heun<double>,
        RK4<double>>;


TYPED_TEST_SUITE_P(RungeKuttaIntegratorTest);


// Simple exponential decay ODE dx/dt = -x with solution x(t) = exp(-t)
Eigen::Matrix<double, 1, 1> decay_function(const Eigen::Matrix<double, 1, 1> &x, double) {
    return -x;
}

// Simple exponential decay ODE dx/dt = -x with solution x(t) = exp(-t)
//auto decay_lambda = [](const Eigen::Matrix<double, 1, 1> &x, double) {
//    return -x;
//};


// Define a system of linear equations
Eigen::Matrix<double, 2, 1> linear_system_function(const Eigen::Matrix<double, 2, 1> &x, double) {
    Eigen::Matrix<double, 2, 1> dx;
    dx[0] = x[1];
    dx[1] = -x[0];
    return dx;
}

//auto linear_system_lambda = [](const Eigen::Matrix<double, 2, 1> &x, double) {
//    Eigen::Matrix<double, 2, 1> dx;
//    dx[0] = x[1];
//    dx[1] = -x[0];
//    return dx;
//};


//TYPED_TEST_P(RungeKuttaIntegratorTest, ExponentialDecayLambda) {
//    this->template TestIntegrator(decay_lambda, std::exp(-t_final), 1);
//}

TYPED_TEST_P(RungeKuttaIntegratorTest, ExponentialDecayFunction) {
    this->template TestIntegrator(decay_function, std::exp(-t_final), 1);
}

//TYPED_TEST_P(RungeKuttaIntegratorTest, LinearSystemLambda) {
//    this->template TestIntegrator(linear_system_lambda, std::cos(t_final), 2);
//}
//
//TYPED_TEST_P(RungeKuttaIntegratorTest, LinearSystemFunction) {
//    this->template TestIntegrator(linear_system_function, std::cos(t_final), 2);
//}

REGISTER_TYPED_TEST_SUITE_P(
        RungeKuttaIntegratorTest,
        ExponentialDecayFunction
        //        ExponentialDecayLambda
);


INSTANTIATE_TYPED_TEST_SUITE_P(IntegratorTypes, RungeKuttaIntegratorTest, IntegratorTypes);


//class EulerMethodTest : public ::testing::Test {
//protected:
//    // Define some common setup, variables here.
//    double h = 0.1;
//    double y0 = 2.0;
//    double u0 = 3.0;
//    double x0 = 0.0;
//
//    virtual void SetUp() {
//        // Setup is called before each test
//    }
//
//    // Define the system of equations based on the given ODE
//    Eigen::VectorXd system_of_equations = [&](const Eigen::VectorXd &x, double t) {
//        Eigen::VectorXd dx(2);
//        dx[0] = x[1];   // y' = u
//        dx[1] = -t * x[1] - x[0]; // u' = -x*u - y
//        return dx;
//    };
//};
//
//TEST_F(EulerMethodTest, ApproximateAtStep) {
//    Eigen::VectorXd initial(2);
//    initial[0] = y0;  // y(0) = 2
//    initial[1] = u0;  // y'(0) = 3
//
//    // Create the Euler integrator for the system of equations
//    RungeKuttaIntegrator<
//            Eigen::VectorXd,
//            decltype(system_of_equations),
//            Euler<double>
//    > eulerIntegrator(system_of_equations);
//
//    // Integrate to get the values at h (0.1)
//    Eigen::VectorXd result = eulerIntegrator.integrate(initial, h, h);
//
//    // Verify the results
//    EXPECT_NEAR(result[0], 2.3, 1e-6);   // y1 = 2.3
//    EXPECT_NEAR(result[1], 2.8, 1e-6);   // u1 = 2.8
//
//    // Integrate again to get the values at 2h (0.2)
//    result = eulerIntegrator.integrate(result, h, h);
//
//    // Verify the results
//    EXPECT_NEAR(result[0], 2.58, 1e-6);   // y2 = 2.58
//    EXPECT_NEAR(result[1], 2.542, 1e-6);  // u2 = 2.542
//}
