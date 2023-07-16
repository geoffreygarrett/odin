#include <Eigen/Dense>
#include <gtest/gtest.h>
#include <odin/vectorize.hpp>
#include <unsupported/Eigen/CXX11/Tensor>


//Eigen::Vector3d acceleration(const Eigen::Vector3d &position) {
//    double sum = position.sum();
//    return {sum, sum, sum};
//}
//
//double potential(const Eigen::Vector3d &position) {
//    return position.sum();
//}
//class PhysicsTest : public ::testing::Test {
//protected:
//    Eigen::Vector3d                          position_vector;
//    Eigen::Matrix<double, 3, Eigen::Dynamic> position_matrix;
//    Eigen::Tensor<double, 3>                 position_tensor_x;
//    Eigen::Tensor<double, 3>                 position_tensor_y;
//    Eigen::Tensor<double, 3>                 position_tensor_z;
//
//    void SetUp() override {
//        position_vector = Eigen::Vector3d(1.0, 2.0, 3.0);
//
//        position_matrix.resize(3, 5);
//        position_matrix << 1, 2, 3, 4, 5,
//                6, 7, 8, 9, 10,
//                11, 12, 13, 14, 15;
//
//        position_tensor_x.resize(5, 5, 5);
//        position_tensor_y.resize(5, 5, 5);
//        position_tensor_z.resize(5, 5, 5);
//        position_tensor_x.setConstant(1.0);
//        position_tensor_y.setConstant(2.0);
//        position_tensor_z.setConstant(3.0);
//    }
//};

TEST(LinspaceTest, BasicTest) {
    // create a linspace from 1 to 5 with 5 elements
    std::vector<double> vec = odin::linspace(1, 5, 5);

    // check size
    ASSERT_EQ(vec.size(), 5);

    // check values
    for (size_t i = 0; i < vec.size(); i++) {
        EXPECT_DOUBLE_EQ(vec[i], i + 1);
    }
}

TEST(LinspaceTest, NegativeRangeTest) {
    // create a linspace from -2 to 2 with 5 elements
    std::vector<double> vec = odin::linspace(-2, 2, 5);

    // check size
    ASSERT_EQ(vec.size(), 5);

    // check values
    for (size_t i = 0; i < vec.size(); i++) {
        EXPECT_DOUBLE_EQ(vec[i], i - 2);
    }
}

TEST(MeshgridTest, BasicTest) {
    // create linspaces
    auto x_vec = odin::eigen::linspace(1., 3., 3);
    auto y_vec = odin::eigen::linspace(4., 6., 3);

    // map to Eigen vectors
    Eigen::Vector<double> x = Eigen::Map<Eigen::Vector<double, Eigen::Dynamic>>(x_vec.data(), x_vec.size());
    Eigen::Vector<double> y = Eigen::Map<Eigen::Vector<double, Eigen::Dynamic>>(y_vec.data(), y_vec.size());

    // create meshgrid
    auto [X, Y] = odin::meshgrid(x, y);

    // check sizes
    ASSERT_EQ(X.dimension(0), 3);
    ASSERT_EQ(X.dimension(1), 3);
    ASSERT_EQ(Y.dimension(0), 3);
    ASSERT_EQ(Y.dimension(1), 3);

    // check values
    for (size_t i = 0; i < x_vec.size(); i++) {
        for (size_t j = 0; j < y_vec.size(); j++) {
            EXPECT_DOUBLE_EQ(X(i, j), x_vec[i]);
            EXPECT_DOUBLE_EQ(Y(i, j), y_vec[j]);
        }
    }
}


Eigen::Vector3d acceleration(const Eigen::Vector3d &position) {
    double sum = position.sum();
    return {sum, sum, sum};
}

double potential(const Eigen::Vector3d &position) {
    return position.sum();
}

Eigen::Vector3d map_index_to_position(const Eigen::Tensor<double, 3> &tensor, int i, int j, int k) {
    return {tensor(i, j, k), tensor(i, j, k), tensor(i, j, k)};
}


TEST(TensorizeTest, TestAcceleration) {
    // Setup data
    Eigen::Tensor<double, 3> input(3, 3, 3);
    input.setRandom();// Fill tensor with random values
    Eigen::Tensor<Eigen::Vector3d, 3> output(3, 3, 3);

    // Wrap function pointer in std::function
    std::function<Eigen::Vector3d(const Eigen::Vector3d &)> accelerationFunc = acceleration;

    // Run tensorize
    odin::tensorize(input, output, accelerationFunc);

    // Check that output matches expected
    for (int i = 0; i < input.dimension(0); ++i) {
        for (int j = 0; j < input.dimension(1); ++j) {
            for (int k = 0; k < input.dimension(2); ++k) {
                Eigen::Vector3d position = map_index_to_position(input, i, j, k);
                EXPECT_EQ(output(i, j, k), acceleration(position));
            }
        }
    }
}

TEST(TensorizeTest, TestPotential) {
    // Setup data
    Eigen::Tensor<double, 3> input(3, 3, 3);
    input.setRandom();// Fill tensor with random values
    Eigen::Tensor<double, 3> output(3, 3, 3);

    // Wrap function pointer in std::function
    std::function<double(const Eigen::Vector3d &)> potentialFunc = potential;

    // Run tensorize
    odin::tensorize(input, output, potentialFunc);

    // Check that output matches expected
    for (int i = 0; i < input.dimension(0); ++i) {
        for (int j = 0; j < input.dimension(1); ++j) {
            for (int k = 0; k < input.dimension(2); ++k) {
                Eigen::Vector3d position = map_index_to_position(input, i, j, k);
                EXPECT_DOUBLE_EQ(output(i, j, k), potential(position));
            }
        }
    }
}

//
//TEST_F(PhysicsTest, TestAccelerationVector) {
//    Eigen::Vector3d expected_acc = Eigen::Vector3d(6.0, 6.0, 6.0);// 1.0 + 2.0 + 3.0
//    Eigen::Vector3d acc          = acceleration(position_vector);
//    EXPECT_TRUE(acc.isApprox(expected_acc));
//}
//
//TEST_F(PhysicsTest, TestPotentialVector) {
//    double expected_pot = 6.0;// 1.0 + 2.0 + 3.0
//    double pot          = potential(position_vector);
//    EXPECT_DOUBLE_EQ(pot, expected_pot);
//}
//
//TEST_F(PhysicsTest, TestAccelerationMatrix) {
//    Eigen::Matrix<double, 3, Eigen::Dynamic> expected_acc(3, 5);
//    expected_acc << 6, 9, 12, 15, 18,// 1+2+3, 2+3+4, 3+4+5, 4+5+6, 5+6+7
//            15, 18, 21, 24, 27,      // 6+7+8, 7+8+9, 8+9+10, 9+10+11, 10+11+12
//            24, 27, 30, 33, 36;      // 11+12+13, 12+13+14, 13+14+15, 14+15+16, 15+16+17
//    auto acc = odin::tensorize(position_matrix, acceleration);
//    EXPECT_TRUE(acc.isApprox(expected_acc));
//}
//
//TEST_F(PhysicsTest, TestPotentialMatrix) {
//    Eigen::VectorXd expected_pot = Eigen::VectorXd(5);
//    expected_pot << 6, 9, 12, 15, 18;// 1+2+3, 2+3+4, 3+4+5, 4+5+6, 5+6+7
//    auto pot = odin::tensorize(position_matrix, potential);
//    EXPECT_TRUE(pot.isApprox(expected_pot));
//}
//
//TEST_F(PhysicsTest, TestAccelerationTensor) {
//    Eigen::Tensor<double, 4> expected_acc(5, 5, 5, 3);
//    expected_acc.setConstant(6.0);// 1.0 + 2.0 + 3.0
//    auto acc = odin::tensorize(position_tensor_x, position_tensor_y, position_tensor_z, acceleration);
//    EXPECT_TRUE(acc.isApprox(expected_acc));
//}
//
//TEST_F(PhysicsTest, TestPotentialTensor) {
//    Eigen::Tensor<double, 3> expected_pot(5, 5, 5);
//    expected_pot.setConstant(6.0);// 1.0 + 2.0 + 3.0
//    auto pot = odin::tensorize(position_tensor_x, position_tensor_y, position_tensor_z, potential);
//    EXPECT_TRUE(pot.isApprox(expected_pot));
//}