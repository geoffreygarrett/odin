#include <include/gtest/gtest.h>
#include <integrators/runge_kutta.hpp>
#include <kalman.hpp>
#include <memory>
#include <random>

class UKFTest : public ::testing::Test {
protected:
    // Earth radius
    static constexpr double R = 6371;

    // Tracking station positions
    static Eigen::Vector2d ts1, ts2, ts3, ts4;

    // Adding some noise
    std::random_device rd;
    std::default_random_engine generator{rd()};
    std::normal_distribution<double> dist{0.0, 0.1};

    // Initial state
//    Eigen::Vector4d true_state{7000, 0, 0, 7.5};
    Eigen::Vector4d true_state{10000, 0, 0, 7.5};

    // Initialize UKF as a pointer
    std::unique_ptr<UnscentedKalmanFilter<Eigen::Vector4d, Eigen::Vector4d, Eigen::Matrix4d, Eigen::Matrix4d>> ukf;

    // Define the functions as member functions
    static Eigen::Vector4d state_transition(const Eigen::Vector4d &x, double dt);

    static Eigen::Vector4d state_dynamics(const Eigen::Vector4d &x);

    static Eigen::Vector4d measurement_model(const Eigen::Vector4d &x);

    static Eigen::Matrix4d process_covariance(const Eigen::Vector4d &, double);

    static Eigen::Matrix4d measurement_covariance(const Eigen::Vector4d &);

    void SetUp() override {
        // Continue with the UKF initialization
        ukf = std::make_unique<UnscentedKalmanFilter<Eigen::Vector4d, Eigen::Vector4d, Eigen::Matrix4d, Eigen::Matrix4d>>(
                UKFTest::state_transition,
                UKFTest::measurement_model,
                UKFTest::process_covariance,
                UKFTest::measurement_covariance,
                0.3,  // alpha
                2.0,  // beta
                1.0   // kappa
        );
        ukf->set_initial_state(true_state);  // no noise added.
        Eigen::Matrix4d initial_covariance = Eigen::Matrix4d::Zero();
//        initial_covariance(2, 2) = 0.001;  // Variance for vx
//        initial_covariance(3, 3) = 0.001;  // Variance for vy
        ukf->set_initial_covariance(initial_covariance);
    }

    void TearDown() override {
        // Perform cleanup actions if needed
    }
};

// Initialize static members
Eigen::Vector2d UKFTest::ts1{-6371, 0};
Eigen::Vector2d UKFTest::ts2{6371, 0};
Eigen::Vector2d UKFTest::ts3{0, 6371};
Eigen::Vector2d UKFTest::ts4{0, -6371};

// Assuming G = 6.67430e-20 km^3 kg^-1 s^-2 (Gravitational constant)
// Assuming m = 5.972e24 kg (Earth's mass)
// Here we're considering the object's mass negligible when compared to Earth's mass

double G = 6.67430e-20;
double m = 5.972e24;

Eigen::Vector4d UKFTest::state_dynamics(const Eigen::Vector4d &x) {
    double r = std::sqrt(std::pow(x[0], 2) + std::pow(x[1], 2));  // distance from origin

    // Calculate accelerations due to gravity
    double a_x = -G * m * x[0] / std::pow(r, 3);
    double a_y = -G * m * x[1] / std::pow(r, 3);

    // Return state derivative
    return Eigen::Vector4d(x[2], x[3], a_x, a_y);
}

Eigen::Vector4d UKFTest::state_transition(const Eigen::Vector4d &x, double dt) {
    Eigen::Vector4d k1 = state_dynamics(x);
    Eigen::Vector4d k2 = state_dynamics(x + 0.5 * dt * k1);
    Eigen::Vector4d k3 = state_dynamics(x + 0.5 * dt * k2);
    Eigen::Vector4d k4 = state_dynamics(x + dt * k3);

    Eigen::Vector4d new_state = x + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4);
    return new_state;
}


Eigen::Vector4d UKFTest::measurement_model(const Eigen::Vector4d &x) {
    Eigen::Vector4d measurement;
    measurement << (Eigen::Vector2d(x[0], x[1]) - ts1).norm(),
            (Eigen::Vector2d(x[0], x[1]) - ts2).norm(),
            (Eigen::Vector2d(x[0], x[1]) - ts3).norm(),
            (Eigen::Vector2d(x[0], x[1]) - ts4).norm();
    return measurement;
}

Eigen::Matrix4d UKFTest::process_covariance(const Eigen::Vector4d &, double) {
    Eigen::Matrix4d Q;
    Q.setZero();
    Q(0, 0) = 0.01;  // Variance for x position in km^2
    Q(1, 1) = 0.01;  // Variance for y position in km^2
    Q(2, 2) = 0.0001; // Variance for x velocity in (km/s)^2
    Q(3, 3) = 0.0001; // Variance for y velocity in (km/s)^2
    return Q;
}

Eigen::Matrix4d UKFTest::measurement_covariance(const Eigen::Vector4d &) {
    return Eigen::Matrix4d::Identity() * 100;
}

#include <fstream>

TEST_F(UKFTest, OrbitEstimation) {
    // Create and open a text file for saving output data
    std::ofstream outputFile("trajectory.csv");

    // Write the headers to the output file
    outputFile << "Iteration,ActualX,ActualY,ActualVx,ActualVy,EstimateX,EstimateY,EstimateVx,EstimateVy,CovXX,CovYY,CovXY\n";

    const double dt = 15;

    // Loop for 1000 steps
    for (int i = 0; i < 10; ++i) {
        // Propagate state with added Gaussian noise
        true_state = state_transition(true_state, dt);

        // Simulate measurement with added Gaussian noise
        Eigen::Vector4d measurement = measurement_model(true_state);
//        measurement[0] += dist(generator);
//        measurement[1] += dist(generator);
//        measurement[2] += dist(generator);
//        measurement[3] += dist(generator);

        // Execute prediction step of UKF
        ukf->predict(dt);

        // Execute update step of UKF with received measurement every 100 steps
        if (i % 400 == 0) {
            ukf->update(measurement, measurement_covariance(measurement));
        }

        // Retrieve estimate state from UKF
        Eigen::Vector4d pred_state = ukf->get_state();
        Eigen::Matrix4d pred_state_covariance = ukf->get_state_covariance();

        // Write the current iteration, actual and estimated state to the output file
        outputFile << i << ","
                   << true_state[0] << "," << true_state[1] << "," << true_state[2] << "," << true_state[3] << ","
                   << pred_state[0] << "," << pred_state[1] << "," << pred_state[2] << "," << pred_state[3] << ","
                   << pred_state_covariance(0,0) << "," << pred_state_covariance(1,1) << "," << pred_state_covariance(0,1)
                   << "\n";
    }

    // Close the output file
    outputFile.close();
}