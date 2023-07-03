#include <Eigen/Dense>
#include <include/gtest/gtest.h>
#include <kalman.hpp>// Adjust this include path as needed

//TEST(SigmaPointCalculatorTest, CalculateSigmaPoints) {
//    int state_dimension = 3;
//    UKFSigmaPointCalculator calc(state_dimension);
//
//    Eigen::VectorXd x(state_dimension);
//    x << 1, 2, 3;
//
//    Eigen::MatrixXd P(state_dimension, state_dimension);
//    P << 1, 0, 0,
//            0, 2, 0,
//            0, 0, 3;
//
//    calc.calculateSigmaPoints(x, P);
//    Eigen::MatrixXd sigma_points = calc.getSigmaPoints();
//
//    // Here we only check if the dimensionality is correct and the first sigma point matches the mean
//    // Further test cases should be implemented based on the specific requirements
//    ASSERT_EQ(sigma_points.rows(), state_dimension);
//    ASSERT_EQ(sigma_points.cols(), 2 * state_dimension + 1);
//    ASSERT_TRUE(sigma_points.col(0).isApprox(x));
//
//    LOG(INFO) << "Sigma points: " << sigma_points;
//}


TEST(SigmaPointCalculatorTest, CalculateSigmaPoints) {
    const int state_dimension = 3;
    UKFSigmaPointCalculator<3> calc;

    Eigen::Vector<double, state_dimension> x;
    x << 1, 2, 3;

    Eigen::Matrix<double, state_dimension, state_dimension> P;
    P << 1, 0, 0,
            0, 2, 0,
            0, 0, 3;

    calc.calculateSigmaPoints(x, P);
    Eigen::MatrixXd sigma_points = calc.get_sigma_points();

    // Here we only check if the dimensionality is correct and the first sigma point matches the mean
    // Further test cases should be implemented based on the specific requirements
    ASSERT_EQ(sigma_points.rows(), state_dimension);
    ASSERT_EQ(sigma_points.cols(), 2 * state_dimension + 1);
    ASSERT_TRUE(sigma_points.col(0).isApprox(x));


    LOG(INFO) << "Sigma points: \n" << sigma_points;
}

class SigmaPointCalculatorFixture : public ::testing::Test {
protected:
    void SetUp() override {
        state_dimension = 3;
        calc = std::make_unique<UKFSigmaPointCalculator<>>(state_dimension);

        x = Eigen::VectorXd(state_dimension);
        x << 1, 2, 3;

        P = Eigen::MatrixXd(state_dimension, state_dimension);
        P << 1, 0, 0,
                0, 2, 0,
                0, 0, 3;

        calc->calculateSigmaPoints(x, P);
    }

    int state_dimension{};
    Eigen::VectorXd x;
    Eigen::MatrixXd P;
    std::unique_ptr<UKFSigmaPointCalculator<>> calc;
};

TEST_F(SigmaPointCalculatorFixture, MeanWeightsSumToOne) {
    // first-order weights must sum to one by definition
    Eigen::VectorXd weights_mean = calc->get_weights_mean();
    double sum_weights_mean = weights_mean.sum();
    ASSERT_NEAR(sum_weights_mean, 1.0, 1e-5);
}

//TEST_F(SigmaPointCalculatorFixture, CovarianceWeightsSumToOne) {
//    // second-order weights must sum to one by definition
//    Eigen::VectorXd weights_covariance = calc->getWeightsCovariance();
//    double sum_weights_covariance = weights_covariance.sum();
//    ASSERT_NEAR(sum_weights_covariance, 1.0, 1e-5);
//}

TEST_F(SigmaPointCalculatorFixture, WeightsArePositive) {
    Eigen::VectorXd weights_mean = calc->get_weights_mean();
    Eigen::VectorXd weights_covariance = calc->get_weights_covariance();

    // Check that all weights are positive
    for (double i: weights_mean) {
        ASSERT_GE(i, 0.0);
    }
    for (double i: weights_covariance) {
        ASSERT_GE(i, 0.0);
    }
}
